# CycFilt
tiny tool for xxx reads filtering



## INSTALL

```shell
git clone https://github.com/comery/CycFilt.git
cd CycFilt
cargo build
# target is : target/debug/cyc_filt
or cargo build --release
# target is : target/release/cyc_filt
```





## Usage

```shell
Usage: cyc_filt [OPTIONS] --input <input_file> --output <output_file> --min-quality <min_quality> --min-length <min_length> --threads <num_threads>

Options:
  -i, --input <input_file>         Input FASTQ file
  -o, --output <output_file>       Output FASTQ file
  -q, --min-quality <min_quality>  Minimum quality score
  -l, --min-length <min_length>    Minimum sequence length
  -t, --threads <num_threads>      Number of threads to use
  -b, --batch-size <batch_size>    Batch size for processing [default: 10000]
  -h, --help                       Print help


e.g.
cyc_filt -i test.fastq.gz -o test.hq.fq.gz -q 7 -l 1000 -t 4
```





## Change logs

- V1.1.1
  - add chunk size for batch reading and processing
- v1.1.0
  - Multi-threads with syn read and write
- V1.0.0
  - multi-threads with async read and write
- v0.1.0
  - single thread for running
- 
