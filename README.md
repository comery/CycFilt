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
Usage: cyc_filt <input_file> <output_file> <min_quality> <min_length>

the input_filt is gzip-format aware, and output fromat is gzip compressed. 

<min_quality> is minimal quality you required

<min_length> is minimal read length you required


e.g.
cyc_filt test.fastq.gz test.hq.fq.gz 7 1000
```





