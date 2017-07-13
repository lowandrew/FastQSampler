# FastQSampler
This project subsamples reads from fastq files to a desired depth.

## Requires:
- Python 3.5 (2.7 should also work - not yet tested).
- Biopython (>=1.67)

## Installation:
* python setup.py build
* python setup.py install

## Use (batch mode):
Use the sampler_wrapper script. Parameters:
* input folder, containing fastq files. Can be either gzipped or plaintext, paired or unpaired. `(-i, --input_folder)`
* output folder, which will hold your subsampled reads. If the folder does not exist, it will be created. `(-o, --output_folder)`
* genome size: Approximate size of the genome you're sampling, in bases. `(-s, --genome_size)`
* coverage depth: Depth of coeverage to sample to. `(-c, --coverage_depth)`

## Optional parameters:
* Gzip output. When included, your output will be gzipped. `(-g, --gzip_output)`
 
