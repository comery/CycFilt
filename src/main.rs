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
    num_cpus: usize,
    batch_size: usize,
    min_quality: f64,
    min_length: usize,
) -> Result<(), IoError> {
    let input_path = std::path::Path::new(input_file);
    let output_path = std::path::Path::new(output_file);

    let input_file = File::open(input_path).expect("Failed to open input file");

    let mut buf = [0; 2];
    File::open(input_path)
        .expect("Failed to open input file again")
        .read_exact(&mut buf)
        .expect("Failed to read first two bytes");
    let reader: Box<dyn BufRead> = if &buf == b"\x1f\x8b" {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(input_file)))
    } else {
        Box::new(BufReader::new(File::open(input_path).expect("Failed to open input file again")))
    };

    let output_file = File::create(output_path).unwrap();
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

        lines.par_chunks(4)
            .enumerate()
            .for_each_with(num_cpus, |_i, chunk| {
                let (_, chunk_ref) = chunk;
                let header = chunk_ref[0].as_ref();
                let sequence = <std::string::String as AsRef<str>>::as_ref(&chunk_ref[1]);
                let quality_value = match get_quality_value(header) {   
                    Ok(val) => val,
                    Err(e) => {
                        eprintln!("Failed to parse quality value from {}: {}", header, e);
                        return;
                    }
                };

                if quality_value >= min_quality && sequence.len() >= min_length {
                    let mut writer_guard = writer.lock().unwrap();
                    for line in chunk_ref {
                        writeln!(&mut *writer_guard, "{}", line).unwrap();
                    }
                } else {
                    let mut filtered_reads_guard = filtered_reads.lock().unwrap();
                    *filtered_reads_guard += 1;
                }

                let mut total_reads_guard = total_reads.lock().unwrap();
                *total_reads_guard += chunk_ref.len();
            });

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


fn get_quality_value(header: &str) -> Result<f64, String> {
    let parts: Vec<&str> = header.split('_').collect();
    if let Some(last_part) = parts.last() {
        last_part.parse::<f64>().map_err(|e| e.to_string())
    } else {
        Err("No underscore found in header".to_string())
    }
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
             .required(true)
             .help("Minimum quality score"))
        .arg(clap::Arg::new("min_length")
             .short('l')
             .long("min-length")
             .required(true)
             .help("Minimum sequence length"))
        .arg(clap::Arg::new("num_cpus")
             .short('c')
             .long("cpus")
             .required(true)
             .help("Number of cpus to use"))
        .arg(clap::Arg::new("batch_size")
             .short('b')
             .long("batch-size")
             .required(false)
             .default_value("10000")
             .help("Batch size for processing"))
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

    let result = filter_fastq_by_quality_and_length(input_file, output_file, num_cpus, batch_size, min_quality, min_length);

    if let Err(e) = result {
        eprintln!("Error: {}", e);
    }
}
