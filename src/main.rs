use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::env;

fn filter_fastq_by_quality_and_length(input: &str, output_file: &str, min_quality: f64, min_length: usize) {
    let input_path = Path::new(input);
    let output_path = Path::new(output_file);

    let input_file = File::open(input_path).expect("Failed to open input file");
    let mut output_file = File::create(output_path).expect("Failed to create output file");

    let reader = BufReader::new(input_file);

    let mut lines = reader.lines();
    let mut total_reads = 0;
    let mut filtered_reads = 0;

    while let Some(Ok(header)) = lines.next() {
        let sequence = lines.next().unwrap().expect("Missing sequence line");
        let plus = lines.next().unwrap().expect("Missing plus line");
        let quality = lines.next().unwrap().expect("Missing quality line");

        total_reads += 1;

        // Debugging output
        // println!("Processing: {}", header);

        // Validate that the header starts with '@ indicating it's a FASTQ header
        if header.starts_with('@') {
            let quality_value = match get_quality_value(&header) {
                Ok(val) => val,
                Err(e) => {
                    eprintln!("Failed to parse quality value from {}: {}", header, e);
                    filtered_reads += 1;
                    continue;
                }
            };

            if quality_value >= min_quality && sequence.len() >= min_length {
                output_file.write_all(header.as_bytes()).unwrap();
                output_file.write_all(b"\n").unwrap();
                output_file.write_all(sequence.as_bytes()).unwrap();
                output_file.write_all(b"\n").unwrap();
                output_file.write_all(plus.as_bytes()).unwrap();
                output_file.write_all(b"\n").unwrap();
                output_file.write_all(quality.as_bytes()).unwrap();
                output_file.write_all(b"\n").unwrap();
            } else {
                filtered_reads += 1;
            }
        } else {
            eprintln!("Invalid FASTQ header: {}", header);
            filtered_reads += 1;
        }
    }

    println!("Total reads processed: {}", total_reads);
    println!("Total reads filtered: {}", filtered_reads);
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
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: {} <input_file> <output_file> <min_quality> <min_length>", args[0]);
        return;
    }

    let input_file = &args[1];
    let output_file = &args[2];
    let min_quality: f64 = args[3].parse().expect("Failed to parse min_quality");
    let min_length: usize = args[4].parse().expect("Failed to parse min_length");
    filter_fastq_by_quality_and_length(input_file, output_file, min_quality, min_length);
}
