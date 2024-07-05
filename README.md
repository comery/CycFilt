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
Usage: cyc_filt --input <FILE> --output <FILE> --quality <FLOAT> --length <INT> --threads <INT>

Options:
  -i, --input <FILE>     Input FASTQ file
  -o, --output <FILE>    Output filtered FASTQ file
  -q, --quality <FLOAT>  Minimum quality threshold
  -l, --length <INT>     Minimum length threshold
  -t, --threads <INT>    Number of threads to use
  -h, --help             Print help
  -V, --version          Print version


e.g.
cyc_filt -i test.fastq.gz -o test.hq.fq.gz -q 7 -l 1000 -t 4
```





