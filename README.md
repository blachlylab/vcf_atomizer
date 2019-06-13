Atomizer
========

Parses VCF files by record into json objects that are written
to stdout. Json objects are newline delimited (jsonl or jsonlines).
Currently the program can parse the ANN field added by SnpEff.
Currently will create a new jsonl object for every unique sample and variant combination.
If ANN field is present in INFO, the atomizer will create a new jsonl object for every 
unique sample, variant, and annotation combination unless the ```-one``` flag is used.

### Requirements
Requires the github.com/brentp/vcfgo package.

### To build
```
make 
```
### To Run
```
./bin/vcf_atomizer sample.vcf.gz
```
or 
```
./bin/vcf_atomizer sample.vcf
```
