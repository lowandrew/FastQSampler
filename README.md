This project subsample reads from fastq files to a desired depth.

Written in Python 3.5

Requires:
Biopython (>=1.67)

Basic Usage:
python3 sampler.py -i <input_fastq> -c <desired_subsample_coverage> -s <size_of_genome>

<input_fastq> is the path to the fastq file you want to subsample. If reads are paired, input both files, separated by a space. These files can be in plain-text or gzipped format.
<desired_subsample_coverage> is the coverage depth you wish to subsample to. If you try to sample more deeply than you can based on your input fastq files, an error will be thrown.
<size_of_genome> is how large you estimate the genome you're sampling to be, in base pairs.

Options:
-o [OUTPUT_FILE] Base name for your output file. Defaults to "out".
-g [GZIP_OUTPUT] Option to gzip your output file(s). True when flag is present, false by default.
