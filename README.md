Atomizer
========

Parses VCF files by record into json objects that are written
to stdout. Json objects are newline delimited (jsonl or jsonlines).
Currently the program can parse the ANN field added by SnpEff.

### Requirements
Requires the github.com/brentp/vcfgo package.

### To build
```
make deps
make all
```
### To Run
```
./bin/atomizer sample.vcf
```
