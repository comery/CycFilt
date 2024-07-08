use std::fs::File;
use std::io::{BufRead, BufReader, Write, Read};
use std::path::Path;
use std::sync::{Arc, Mutex, Barrier};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use clap::{Arg, Command};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use crossbeam::channel::bounded;
use std::io::Error as IoError;

fn filter_fastq_by_quality_and_length(
    input_file: &str,
    output_file: &str,
    min_quality: f64,
    min_length: usize,
    num_threads: usize,
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

    let output_file = File::create(output_path).expect("Failed to create output file");
    let writer = GzEncoder::new(output_file, Compression::default());

    let writer = Arc::new(Mutex::new(writer));
    let barrier = Arc::new(Barrier::new(num_threads + 1));

    let (sender, receiver) = bounded::<(String, String, String, String)>(num_threads);

    let total_reads = Arc::new(AtomicUsize::new(0));
    let filtered_reads = Arc::new(AtomicUsize::new(0));

    for _ in 0..num_threads {
        let receiver = receiver.clone();
        let writer = Arc::clone(&writer);
        let barrier = Arc::clone(&barrier);
        let total_reads = Arc::clone(&total_reads);
        let filtered_reads = Arc::clone(&filtered_reads);

        thread::spawn(move || {
            while let Ok((header, sequence, plus, quality)) = receiver.recv() {
                total_reads.fetch_add(1, Ordering::Relaxed);

                let quality_value = match get_quality_value(&header) {
                    Ok(val) => val,
                    Err(e) => {
                        eprintln!("Failed to parse quality value from {}: {}", header, e);
                        filtered_reads.fetch_add(1, Ordering::Relaxed);
                        continue;
                    }
                };

                if quality_value >= min_quality && sequence.len() >= min_length {
                    let mut writer = writer.lock().unwrap();
                    writer.write_all(header.as_bytes()).unwrap();
                    writer.write_all(b"\n").unwrap();
                    writer.write_all(sequence.as_bytes()).unwrap();
                    writer.write_all(b"\n").unwrap();
                    writer.write_all(plus.as_bytes()).unwrap();
                    writer.write_all(b"\n").unwrap();
                    writer.write_all(quality.as_bytes()).unwrap();
                    writer.write_all(b"\n").unwrap();
                } else {
                    filtered_reads.fetch_add(1, Ordering::Relaxed);
                }
            }
            barrier.wait();
        });
    }

    let mut lines = reader.lines();

    while let Some(Ok(header)) = lines.next() {
        let sequence = lines.next().unwrap_or(Ok(String::new())).expect("Missing sequence line");
        let plus = lines.next().unwrap_or(Ok(String::new())).expect("Missing plus line");
        let quality = lines.next().unwrap_or(Ok(String::new())).expect("Missing quality line");


        // println!("Processing: {}", header);

        if header.starts_with('@') {
            sender.send((header, sequence, plus, quality)).unwrap();
        } else {
            eprintln!("Invalid FASTQ header: {}", header);
            filtered_reads.fetch_add(1, Ordering::Relaxed);
        }
    }

    drop(sender);
    barrier.wait();

    // Instead of unwrapping, we'll just lock the mutex and finish the compression
    {
        let mut writer = writer.lock().map_err(|_| IoError::new(std::io::ErrorKind::Other, "Failed to lock writer"))?;
        writer.try_finish()?;
    }

    let total_reads = total_reads.load(Ordering::Relaxed);
    let filtered_reads = filtered_reads.load(Ordering::Relaxed);

    println!("Total reads processed: {}", total_reads);
    println!("Total reads filtered: {}", filtered_reads);

    Ok(())
}

fn get_quality_value(header: &str) -> Result<f64, String> {
    let parts: Vec<&str> = header.split('_').collect();
    if let Some(last_part) = parts.last() {
        last_part.parse::<f64>().map_err(|e: std::num::ParseFloatError| e.to_string())
    } else {
        Err("No underscore found in header".to_string())
    }
}

fn main() {
    let matches = Command::new("fastq_filter")
        .version("1.0")
        .about("Filters FASTQ files based on quality and length")
        .arg(Arg::new("input_file")
            .short('i')
            .long("input")
            .value_name("FILE")
            .required(true)
            .help("Input FASTQ file"))
        .arg(Arg::new("output_file")
            .short('o')
            .long("output")
            .value_name("FILE")
            .required(true)
            .help("Output filtered FASTQ file"))
        .arg(Arg::new("min_quality")
            .short('q')
            .long("quality")
            .value_name("FLOAT")
            .required(true)
            .help("Minimum quality threshold"))
        .arg(Arg::new("min_length")
            .short('l')
            .long("length")
            .value_name("INT")
            .required(true)
            .help("Minimum length threshold"))
        .arg(Arg::new("num_threads")
            .short('t')
            .long("threads")
            .value_name("INT")
            .required(true)
            .help("Number of threads to use"))
        .get_matches();

    let input_file = matches.get_one::<String>("input_file").unwrap();
    let output_file = matches.get_one::<String>("output_file").unwrap();
    let min_quality: f64 = matches.get_one::<String>("min_quality").unwrap().parse().expect("Failed to parse min_quality");
    let min_length: usize = matches.get_one::<String>("min_length").unwrap().parse().expect("Failed to parse min_length");
    let num_threads: usize = matches.get_one::<String>("num_threads").unwrap().parse().expect("Failed to parse num_threads");

    match filter_fastq_by_quality_and_length(input_file, output_file, min_quality, min_length, num_threads) {
        Ok(_) => println!("Filtering completed successfully."),
        Err(e) => eprintln!("An error occurred: {}", e),
    }
}
