use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use tokio::sync::mpsc;
use clap::{Arg, Command};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use crossbeam::channel;

// 异步写操作函数
async fn async_write(
    writer_arc: Arc<Mutex<GzEncoder<File>>>,
    mut receiver: mpsc::Receiver<(String, String, String, String)>
) {
    while let Some((header, sequence, plus, quality)) = receiver.recv().await {
        let mut writer = writer_arc.lock().unwrap();
        writer.write_all(header.as_bytes()).unwrap();
        writer.write_all(sequence.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        writer.write_all(plus.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
        writer.write_all(quality.as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();
    }

    // Move writer.finish() to be called outside the Mutex guard
    let writer = Arc::try_unwrap(writer_arc)
        .expect("Failed to unwrap Arc")
        .into_inner()
        .expect("Failed to unwrap Mutex");
    writer.finish().unwrap();
}

fn filter_fastq_by_quality_and_length(
    input_file: &str,
    output_file: &str,
    min_quality: f64,
    min_length: usize,
    num_threads: usize,
) {
    let input_path = Path::new(input_file);
    let output_path = Path::new(output_file);

    let mut input_file = File::open(input_path).expect("Failed to open input file");
    let mut buf = [0; 2];
    input_file.read_exact(&mut buf).expect("Failed to read first two bytes");

    // Input needs to be reopened because first 2 bytes were read to identify if gzip
    let input_file = File::open(input_path).expect("Failed to reopen input file");
    let reader: Box<dyn BufRead> = if &buf == b"\x1f\x8b" {
        Box::new(BufReader::with_capacity(16 * 1024, GzDecoder::new(input_file)))  // 增大缓冲区
    } else {
        Box::new(BufReader::with_capacity(16 * 1024, input_file))
    };

    let output_file = File::create(output_path).expect("Failed to create output file");
    let writer = GzEncoder::new(output_file, Compression::default());

    let writer_arc = Arc::new(Mutex::new(writer));

    let (sender, receiver) = channel::bounded::<(String, String, String, String)>(num_threads *2);  // 增加通道容量

    let (async_sender, async_receiver) = mpsc::channel(100);

    let writer_arc_clone = Arc::clone(&writer_arc);

    tokio::spawn(async move {
        async_write(writer_arc_clone, async_receiver).await;
    });

    let handles: Vec<_> = (0..num_threads).map(|_| {
        let receiver = receiver.clone();
        let async_sender = async_sender.clone();

        std::thread::spawn(move || {
            while let Ok((header, sequence, plus, quality)) = receiver.recv() {
                let quality_value = match get_quality_value(&header) {
                    Ok(val) => val,
                    Err(e) => {
                        eprintln!("Failed to parse quality value from {}: {}", header, e);
                        continue;
                    }
                };

                if quality_value >= min_quality && sequence.len() >= min_length {
                    async_sender.blocking_send((header, sequence, plus, quality)).unwrap();
                }
            }
        })
    }).collect();

    let mut total_reads = 0;
    let mut filtered_reads = 0;

    let mut lines = reader.lines();

    while let Some(Ok(header)) = lines.next() {
        let sequence = lines.next().unwrap_or(Ok(String::new())).expect("Missing sequence line");
        let plus = lines.next().unwrap_or(Ok(String::new())).expect("Missing plus line");
        let quality = lines.next().unwrap_or(Ok(String::new())).expect("Missing quality line");

        total_reads += 1;

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
            filtered_reads += 1;
            sender.send((header, sequence, plus, quality)).unwrap();
        }
    }

    drop(sender);  // Close the sending channel to finish threads

    for handle in handles {
        handle.join().unwrap();
    }

    println!("Total reads processed: {}", total_reads);
    println!("Total reads passing filters: {}", filtered_reads);
}

fn get_quality_value(header: &str) -> Result<f64, String> {
    let parts: Vec<&str> = header.split('_').collect();
    if let Some(last_part) = parts.last() {
        last_part.parse::<f64>().map_err(|e| e.to_string())
    } else {
        Err("No underscore found in header".to_string())
    }
}

#[tokio::main]
async fn main() {
    let matches = Command::new("fastq_filter")
        .version("1.0")
        .about("Filters FASTQ files based on quality and length")
        .arg(Arg::new("input_file")
            .short('i')
            .long("input")
            .required(true)
            .value_name("FILE")
            .help("Input FASTQ file"))
        .arg(Arg::new("output_file")
            .short('o')
            .long("output")
            .required(true)
            .value_name("FILE")
            .help("Output filtered FASTQ file"))
        .arg(Arg::new("min_quality")
            .short('q')
            .long("quality")
            .required(true)
            .value_name("FLOAT")
            .help("Minimum quality threshold"))
        .arg(Arg::new("min_length")
            .short('l')
            .long("length")
            .required(true)
            .value_name("INT")
            .help("Minimum length threshold"))
        .arg(Arg::new("num_threads")
            .short('t')
            .long("threads")
            .required(true)
            .value_name("INT")
            .help("Number of threads to use"))
        .get_matches();

    let input_file = matches.get_one::<String>("input_file").unwrap();
    let output_file = matches.get_one::<String>("output_file").unwrap();
    let min_quality: f64 = matches.get_one::<String>("min_quality").unwrap().parse().expect("Failed to parse min_quality");
    let min_length: usize = matches.get_one::<String>("min_length").unwrap().parse().expect("Failed to parse min_length");
    let num_threads: usize = matches.get_one::<String>("num_threads").unwrap().parse().expect("Failed to parse num_threads");

    filter_fastq_by_quality_and_length(input_file, output_file, min_quality, min_length, num_threads);
}
