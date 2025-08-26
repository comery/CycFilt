use std::sync::{Arc, Mutex};
use std::io::{BufReader, BufRead, Read, BufWriter, Write};
use std::fs::File;
// use std::path::Path;
// use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Error as IoError;
use rayon::prelude::*;
use num_cpus;

fn filter_fastq_by_quality_and_length(
    input_file: &str,
    output_file: &str,
    _num_cpus: usize,
    batch_size: usize,
    min_quality: f64,
    min_length: usize,
    adapter_sequence: Option<&str>,
    min_adapter_match: usize,
    max_mismatches: usize,
    max_indels: usize,
    debug_mode: bool,
) -> Result<(), IoError> {
    let input_path = std::path::Path::new(input_file);
    let output_path = std::path::Path::new(output_file);

    let mut buf = [0; 2];
    let mut input_file_for_check = File::open(input_path).expect("Failed to open input file");
    input_file_for_check
        .read_exact(&mut buf)
        .expect("Failed to read first two bytes");
    let reader: Box<dyn BufRead> = if &buf == b"\x1f\x8b" {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(File::open(input_path)?)))
    } else {
        Box::new(BufReader::new(File::open(input_path)?))
    };

    let output_file = File::create(output_path)?;
    let writer = Arc::new(Mutex::new(GzEncoder::new(BufWriter::new(output_file), Compression::default())));

    let total_reads = Arc::new(Mutex::new(0));
    let filtered_reads = Arc::new(Mutex::new(0));

    let mut batch_start = 0;

    let mut lines_iter = reader.lines();
    loop {
        let lines: Vec<_> = lines_iter.by_ref().take(batch_size * 4).collect::<Result<Vec<_>, _>>()?;
        if lines.is_empty() {
            break;
        }

        let batch_end = batch_start + lines.len() / 4;

        let writer_clone = Arc::clone(&writer);
        let (local_total, local_filtered, output_lines) = lines.par_chunks(4)
            .fold(
                || (0, 0, Vec::new()),
                |(total, mut filtered, mut output_lines), chunk| {
                    let header = &chunk[0];
                    let sequence = &chunk[1];
                    let quality_line = &chunk[3];
                    let quality_value = match get_quality_value(header) {   
                        Ok(val) => val,
                        Err(e) => {
                            if debug_mode {
                                eprintln!("DEBUG: Failed to parse quality value from {}: {}", header, e);
                            }
                            return (total + chunk.len() / 4, filtered + 1, output_lines);
                        }
                    };

                    if quality_value >= min_quality && sequence.len() >= min_length {
                        if let Some(adapter) = adapter_sequence {
                            let processed_seqs = process_adapter_sequence(
                                header, sequence, quality_line, adapter, 
                                min_adapter_match, max_mismatches, max_indels, debug_mode
                            );
                            
                            for (processed_header, processed_seq, processed_qual) in processed_seqs {
                                if processed_seq.len() >= min_length {
                                    output_lines.push(processed_header);
                                    output_lines.push(processed_seq);
                                    output_lines.push("+".to_string());
                                    output_lines.push(processed_qual);
                                } else {
                                    if debug_mode {
                                        eprintln!("DEBUG: Filtered {} - trimmed length {} < {}", processed_header, processed_seq.len(), min_length);
                                    }
                                    filtered += 1;
                                }
                            }
                        } else {
                            output_lines.extend(chunk.iter().map(|s| s.clone()));
                        }
                    } else {
                        if debug_mode {
                            if quality_value < min_quality {
                                eprintln!("DEBUG: Filtered {} - quality {} < {}", header, quality_value, min_quality);
                            } else {
                                eprintln!("DEBUG: Filtered {} - length {} < {}", header, sequence.len(), min_length);
                            }
                        }
                        filtered += 1;
                    }
                    
                    (total + chunk.len() / 4, filtered, output_lines)
                }
            )
            .reduce(
                || (0, 0, Vec::new()),
                |(total1, filtered1, mut lines1), (total2, filtered2, lines2)| {
                    lines1.extend(lines2);
                    (total1 + total2, filtered1 + filtered2, lines1)
                }
            );

        if !output_lines.is_empty() {
            let mut writer_guard = writer_clone.lock().unwrap();
            for line in &output_lines {
                writeln!(&mut *writer_guard, "{}", line)?;
            }
        }

        {
            let mut total_reads_guard = total_reads.lock().unwrap();
            *total_reads_guard += local_total;
        }
        {
            let mut filtered_reads_guard = filtered_reads.lock().unwrap();
            *filtered_reads_guard += local_filtered;
        }

        batch_start = batch_end;
    }

    let mut writer_guard = writer.lock().unwrap();
    writer_guard.flush()?;

    let total_reads = *total_reads.lock().unwrap();
    let filtered_reads = *filtered_reads.lock().unwrap();
    println!("Total reads: {}", total_reads);
    println!("Filtered reads: {}", filtered_reads);

    Ok(())
}


fn get_quality_value(header: &String) -> Result<f64, String> {
    let parts: Vec<&str> = header.split('_').collect();
    if let Some(last_part) = parts.last() {
        // Handle cases where there might be additional text after the quality value
        let quality_str = last_part.split_whitespace().next().unwrap_or(last_part);
        quality_str.parse::<f64>().map_err(|e| e.to_string())
    } else {
        Err("No underscore found in header".to_string())
    }
}

fn smith_waterman_align(
    sequence: &str,
    adapter: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    max_mismatches: usize,
    max_indels: usize,
) -> Option<(usize, usize, usize, usize)> {
    let seq_chars: Vec<char> = sequence.chars().collect();
    let adapter_chars: Vec<char> = adapter.chars().collect();
    let m = seq_chars.len();
    let n = adapter_chars.len();
    
    if m == 0 || n == 0 {
        return None;
    }
    
    let mut matrix = vec![vec![0; n + 1]; m + 1];
    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;
    
    for i in 1..=m {
        for j in 1..=n {
            let match_val = if seq_chars[i-1] == adapter_chars[j-1] {
                matrix[i-1][j-1] + match_score
            } else {
                matrix[i-1][j-1] + mismatch_penalty
            };
            
            let delete = matrix[i-1][j] + gap_penalty;
            let insert = matrix[i][j-1] + gap_penalty;
            
            matrix[i][j] = 0.max(match_val).max(delete).max(insert);
            
            if matrix[i][j] > max_score {
                max_score = matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    if max_score <= 0 {
        return None;
    }
    
    let mut i = max_i;
    let mut j = max_j;
    let mut mismatches = 0;
    let mut indels = 0;
    
    while i > 0 && j > 0 && matrix[i][j] > 0 {
        if matrix[i][j] == matrix[i-1][j-1] + 
            if seq_chars[i-1] == adapter_chars[j-1] { match_score } else { mismatch_penalty } {
            if seq_chars[i-1] != adapter_chars[j-1] {
                mismatches += 1;
            }
            i -= 1;
            j -= 1;
        } else if matrix[i][j] == matrix[i-1][j] + gap_penalty {
            indels += 1;
            i -= 1;
        } else if matrix[i][j] == matrix[i][j-1] + gap_penalty {
            indels += 1;
            j -= 1;
        } else {
            break;
        }
    }
    
    if mismatches <= max_mismatches && indels <= max_indels {
        Some((i, max_i, j, max_j))
    } else {
        None
    }
}

fn detect_adapter_position(
    sequence: &str,
    adapter: &str,
    min_match: usize,
    max_mismatches: usize,
    max_indels: usize,
) -> Option<usize> {
    if adapter.len() < min_match {
        return None;
    }
    
    if let Some((start_i, _end_i, start_j, end_j)) = smith_waterman_align(
        sequence, adapter, 2, -1, -2, max_mismatches, max_indels
    ) {
        let aligned_length = end_j - start_j;
        if aligned_length >= min_match {
            return Some(start_i);
        }
    }
    
    None
}

fn process_adapter_sequence(
    header: &str,
    sequence: &str,
    quality: &str,
    adapter: &str,
    min_match: usize,
    max_mismatches: usize,
    max_indels: usize,
    debug_mode: bool,
) -> Vec<(String, String, String)> {
    let mut results = Vec::new();
    
    if let Some(pos) = detect_adapter_position(sequence, adapter, min_match, max_mismatches, max_indels) {
        if debug_mode {
            eprintln!("DEBUG: Adapter found in {} at position {}", header, pos);
        }
        
        // Always split at adapter position
        let part1_seq = &sequence[..pos];
        let part1_qual = &quality[..pos];
        if !part1_seq.is_empty() {
            results.push((
                format!("{}_part1", header),
                part1_seq.to_string(),
                part1_qual.to_string()
            ));
        }
        
        let part2_seq = &sequence[pos..];
        let part2_qual = &quality[pos..];
        if !part2_seq.is_empty() {
            results.push((
                format!("{}_part2", header),
                part2_seq.to_string(),
                part2_qual.to_string()
            ));
        }
    } else {
        // No adapter found, keep original
        results.push((
            header.to_string(),
            sequence.to_string(),
            quality.to_string()
        ));
    }
    
    results
}


fn main() {
    // let default_batch_size: usize = 10000;
    let matches = clap::Command::new("fastq-filter")
        .arg(clap::Arg::new("input_file")
             .short('i')
             .long("input")
             .required(true)
             .help("Input FASTQ file"))
        .arg(clap::Arg::new("output_file")
             .short('o')
             .long("output")
             .required(true)
             .help("Output FASTQ file"))
        .arg(clap::Arg::new("min_quality")
             .short('q')
             .long("min-quality")
             .required(false)
             .default_value("7.0")
             .help("Minimum quality score"))
        .arg(clap::Arg::new("min_length")
             .short('l')
             .long("min-length")
             .required(false)
             .default_value("1000")
             .help("Minimum sequence length"))
        .arg(clap::Arg::new("num_cpus")
             .short('c')
             .long("cpus")
             .required(false)
             .default_value("2")
             .help("Number of cpus to use"))
        .arg(clap::Arg::new("batch_size")
             .short('b')
             .long("batch-size")
             .required(false)
             .default_value("10000")
             .help("Batch size for processing"))
        .arg(clap::Arg::new("adapter")
             .short('a')
             .long("adapter")
             .required(false)
             .help("Adapter sequence to detect and remove"))
        .arg(clap::Arg::new("min_adapter_match")
             .short('m')
             .long("min-adapter-match")
             .required(false)
             .default_value("10")
             .help("Minimum adapter match length"))
        .arg(clap::Arg::new("max_mismatches")
             .short('x')
             .long("max-mismatches")
             .required(false)
             .default_value("2")
             .help("Maximum allowed mismatches in adapter alignment"))
        .arg(clap::Arg::new("max_indels")
             .short('d')
             .long("max-indels")
             .required(false)
             .default_value("1")
             .help("Maximum allowed indels in adapter alignment"))
        .arg(clap::Arg::new("debug")
             .short('D')
             .long("debug")
             .required(false)
             .action(clap::ArgAction::SetTrue)
             .help("Enable debug output with detailed filtering information"))
        .get_matches();

    let input_file = matches.get_one::<String>("input_file").unwrap();
    let output_file = matches.get_one::<String>("output_file").unwrap();

    // Check if the input file can be opened
    let input_file_path = std::path::Path::new(input_file);
    if !input_file_path.exists() {
        eprintln!("Error: input file '{}' does not exist.", input_file);
        std::process::exit(1);
    }

    if !input_file_path.is_file() {
        eprintln!("Error: '{}' is not a file.", input_file);
        std::process::exit(1);
    }

    let min_quality_str = matches.get_one::<String>("min_quality").unwrap();
    let min_quality: f64 = match min_quality_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'min_quality'. Expected a floating-point number.");
            std::process::exit(1);
        }
    };

    let min_length_str = matches.get_one::<String>("min_length").unwrap();
    let min_length: usize = match min_length_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'min_length'. Expected a positive integer.");
            std::process::exit(1);
        }
    };

    let all_cpus = num_cpus::get(); // Get the number of available CPUs
    let num_cpus_str = matches.get_one::<String>("num_cpus").unwrap();
    let num_cpus: usize = match num_cpus_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'num_cpus'. Expected a positive integer.");
            std::process::exit(1);
        }
    };
    if num_cpus > all_cpus {
        eprintln!("Error: inalid value for 'num_cpus' beacuse it exceed to all avaliable cpus {}", all_cpus);
        std::process::exit(1);
    }

    let batch_size_str = matches.get_one::<String>("batch_size").unwrap();
    let batch_size: usize = match batch_size_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'batch_size'. Expected a positive integer.");
            std::process::exit(1);
        }
    };

    let adapter_sequence = matches.get_one::<String>("adapter").map(|s| s.as_str());
    
    let min_adapter_match_str = matches.get_one::<String>("min_adapter_match").unwrap();
    let min_adapter_match: usize = match min_adapter_match_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'min_adapter_match'. Expected a positive integer.");
            std::process::exit(1);
        }
    };


    let max_mismatches_str = matches.get_one::<String>("max_mismatches").unwrap();
    let max_mismatches: usize = match max_mismatches_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'max_mismatches'. Expected a positive integer.");
            std::process::exit(1);
        }
    };

    let max_indels_str = matches.get_one::<String>("max_indels").unwrap();
    let max_indels: usize = match max_indels_str.parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error: invalid value for 'max_indels'. Expected a positive integer.");
            std::process::exit(1);
        }
    };

    let debug_mode = matches.get_flag("debug");

    let result = filter_fastq_by_quality_and_length(
        input_file, 
        output_file, 
        num_cpus, 
        batch_size, 
        min_quality, 
        min_length,
        adapter_sequence,
        min_adapter_match,
        max_mismatches,
        max_indels,
        debug_mode
    );

    if let Err(e) = result {
        eprintln!("Error: {}", e);
    }
}
