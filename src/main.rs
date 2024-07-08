use std::sync::{Arc, Mutex};
use std::io::{BufReader, BufRead, Read, BufWriter, Write};
use std::fs::File;
use std::path::Path;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Error as IoError;

fn filter_fastq_by_quality_and_length(
    input_file: &str,
    output_file: &str,
    num_threads: usize,
    batch_size: usize,
    min_quality: f64,
    min_length: usize,
) -> Result<(), IoError> {

    let input_path = Path::new(input_file);
    let output_path = Path::new(output_file);

    let input_file = File::open(input_path).expect("Failed to open input file");

    let mut buf = [0; 2];
    File::open(input_path)
        .expect("Failed to open input file again")
        .read_exact(&mut buf)
        .expect("Failed to read first two bytes");
    let reader: Box<dyn BufRead> = if &buf == b"\x1f\x8b" {
        Box::new(BufReader::new(GzDecoder::new(input_file)))
    } else {
        Box::new(BufReader::new(File::open(input_path).expect("Failed to open input file again")))
    };

    // Create the GzEncoder once
    let output_file = File::create(output_path).unwrap();
    let writer = Arc::new(Mutex::new(GzEncoder::new(BufWriter::new(output_file), Compression::default())));

    let total_reads = Arc::new(Mutex::new(0));
    let filtered_reads = Arc::new(Mutex::new(0));

    let mut batch_start = 0;

    // Iterate directly over the reader's lines
    let mut lines_iter = reader.lines(); 
    //while let Ok(lines) = lines_iter.by_ref().take(batch_size * 4).collect::<Result<Vec<_>, _>>() {
    loop { // Use loop instead of while let
        let lines: Vec<_> = lines_iter.by_ref().take(batch_size * 4).collect::<Result<Vec<_>, _>>()?; // Collect lines and handle potential errors
        if lines.is_empty() { // Break if no lines were read
            break;
        }

        let batch_end = batch_start + lines.len() / 4;

        let chunk_size = lines.len() / (4 * num_threads); // Correct chunk size calculation

        let threads: Vec<_> = (0..num_threads)
            .map(|thread_id| {
                let start = thread_id * chunk_size * 4; // Correct starting index
                let end = std::cmp::min((thread_id + 1) * chunk_size * 4, lines.len()); // Correct ending index
                let lines_clone = lines.clone();
                let total_reads_arc = Arc::clone(&total_reads);
                let filtered_reads_arc = Arc::clone(&filtered_reads);
                let writer_arc = Arc::clone(&writer);

                std::thread::spawn(move || {
                    let mut thread_total_reads = 0;
                    let mut thread_filtered_reads = 0;

                    for i in start..end {
                        if (i % 4) == 0 { // Only process headers
                            let header = <std::string::String as AsRef<str>>::as_ref(&lines_clone[i]);
                            let sequence = <std::string::String as AsRef<str>>::as_ref(&lines_clone[i + 1]);
                            // println!("Processing: {}", header);
                            // Check if the sequence passes the quality and length filter
                            let quality_value = match get_quality_value(&header) {
                                Ok(val) => val,
                                Err(e) => {
                                    eprintln!("Failed to parse quality value from {}: {}", header, e);
                                    continue;
                                }
                            };

                            if quality_value >= min_quality && sequence.len() >= min_length {
                                // Acquire lock on the writer
                                let mut writer_guard = writer_arc.lock().unwrap();
                                for j in 0..4 {
                                    writeln!(&mut *writer_guard, "{}", <std::string::String as AsRef<str>>::as_ref(&lines_clone[i + j])).unwrap();
                                }
                            } else {
                                thread_filtered_reads += 1;
                            }
                            thread_total_reads += 1;
                        }
                    }
                    // println!("Processing: {}", batch_end);
                    let mut total_reads_guard = total_reads_arc.lock().unwrap();
                    *total_reads_guard += thread_total_reads;

                    let mut filtered_reads_guard = filtered_reads_arc.lock().unwrap();
                    *filtered_reads_guard += thread_filtered_reads;
                })
            })
            .collect();

        for thread in threads {
            thread.join().unwrap();
        }

        batch_start = batch_end;
    }
     // Flush the GzEncoder to ensure all data is written to the file
    let mut writer_guard = writer.lock().unwrap();
    writer_guard.flush()?;

    let total_reads = *total_reads.lock().unwrap();
    let filtered_reads = *filtered_reads.lock().unwrap();
    println!("Total reads: {}", total_reads);
    println!("Filtered reads: {}", filtered_reads);

    Ok(())

}

fn get_quality_value(header: &str) -> Result<f64, String> {
    let parts: Vec<&str> = header.split('_').collect();
    if let Some(last_part) = parts.last() {
        last_part.parse::<f64>().map_err(|e| e.to_string())
    } else {
        Err("No underscore found in header".to_string())
    }
}


fn main() {
    let default_batch_size: usize = 10000;
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
             .required(true)
             .help("Minimum quality score"))
        .arg(clap::Arg::new("min_length")
             .short('l')
             .long("min-length")
             .required(true)
             .help("Minimum sequence length"))
        .arg(clap::Arg::new("num_threads")
             .short('t')
             .long("threads")
             .required(true)
             .help("Number of threads to use"))
        .arg(clap::Arg::new("batch_size")
             .short('b')
             .long("batch-size")
             .required(false)
             .default_value(&default_batch_size.to_string())
             .help("Batch size for processing"))
        .get_matches();

    let input_file = matches.get_one::<String>("input_file").unwrap();
    let output_file = matches.get_one::<String>("output_file").unwrap();
    let min_quality: f64 = matches.get_one::<String>("min_quality").unwrap().parse().unwrap();
    let min_length: usize = matches.get_one::<String>("min_length").unwrap().parse().unwrap();
    let num_threads: usize = matches.get_one::<String>("num_threads").unwrap().parse().unwrap();
    let batch_size: usize = matches.get_one::<String>("batch_size").unwrap().parse().unwrap();

    let result = filter_fastq_by_quality_and_length(input_file, output_file, num_threads, batch_size, min_quality, min_length);

    if let Err(e) = result {
        eprintln!("Error: {}", e);
    }
}
