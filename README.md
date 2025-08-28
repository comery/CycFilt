# CycFilt
tiny tool for CycloneSEQ reads filtering


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
Usage: cyc_filt [OPTIONS] --input <input_file> --output <output_file>

Options:
  -i, --input <input_file>
          Input FASTQ file
  -o, --output <output_file>
          Output FASTQ file
  -q, --min-quality <min_quality>
          Minimum quality score [default: 7.0]
  -l, --min-length <min_length>
          Minimum sequence length [default: 1000]
  -c, --cpus <num_cpus>
          Number of cpus to use [default: 2]
  -b, --batch-size <batch_size>
          Batch size for processing [default: 10000]
  -a, --adapter <adapter>
          Adapter sequence to detect and remove
  -m, --min-adapter-match <min_adapter_match>
          Minimum adapter match length [default: 10]
  -x, --max-mismatches <max_mismatches>
          Maximum allowed mismatches in adapter alignment [default: 2]
  -d, --max-indels <max_indels>
          Maximum allowed indels in adapter alignment [default: 1]
  -D, --debug
          Enable debug output with detailed filtering information
  -h, --help
          Print help


e.g.
cyc_filt -i test.fastq.gz -o test.hq.fq.gz -q 7 -l 1000 -t 4
cyc_filt -i test.fastq.gz -o test.hq.fq.gz -a GGGTGACAGAGCAAGACCCTGTCTCAGAA  -x 3 -d 1  -D

```





## Change logs

- V2.0.1 # 20250829
   - improve the adapter sequence identification to handle reverse complement cases

- V2.0.0 # 20250826
   - Improve performance
   - Add the function of adapter detecting and trimming

- V1.2.0, # 20240708

  - Multi CPUs instead of threads

- V1.1.1

  - add chunk size for batch reading and processing

- v1.1.0

  - Multi-threads with syn read and write

- V1.0.0

  - multi-threads with async read and write

- v0.1.0

  - single thread for running

  
