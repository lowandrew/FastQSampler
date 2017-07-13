#!/usr/bin/env python

import os


def parse_fastq_folder(fastqfolder):
    import glob
    print("Parsing fastq directory...")
    fastq_pairs = list()
    fastq_singles = list()
    # Get a list of all fastq files. For some reason, having/not having the slash doesn't seem to matter on the
    # fastqfolder argument. These should be all the common extensions
    fastq_files = glob.glob(fastqfolder + "/*.fastq*")
    fastq_files += glob.glob(fastqfolder + "/*.fq*")
    for name in fastq_files:
        # If forward and reverse reads are present, put them in a list of paired files.
        # May need to add support for other naming conventions too. Supports both _R1 and _1 type conventions.
        if "R1" in name and os.path.isfile(name.replace("R1", "R2")):
            fastq_pairs.append([name, name.replace("R1", "R2")])
        # If not paired, put in a list of single fastq files.
        elif "_1" in name and os.path.isfile(name.replace("_1", "_2")):
            fastq_pairs.append([name, name.replace("_1", "_2")])
        # This should probably get changed to an else.
        elif "R1" not in name and "R2" not in name and "_1" not in name and "_2" not in name:
            fastq_singles.append(name)

    return fastq_pairs, fastq_singles


def create_output_folder(outdir):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", help="Folder with reads to subsample in FASTQ format.", required=True)
    parser.add_argument("-c", "--coverage_depth", type=float, required=True, help="Coverage depth to subsample to.")
    parser.add_argument("-s", "--genome_size", type=int, required=True, help="Expected genome size in bp.")
    parser.add_argument("-g", "--gzip_output", default=False, action="store_true", help="When included"
                                                                                        ", gzips output.")
    parser.add_argument("-o", "--output_folder", required=True, help="Output folder to store your subsampled reads.")
    arguments = parser.parse_args()
    create_output_folder(arguments.output_folder)
    pairs, singles = parse_fastq_folder(arguments.input_folder)
    for pair in pairs:
        outname = pair[0].split("/")[-1]
        outname = outname.split("_")[0] + "_subsampled" + str(arguments.coverage_depth)
        outname = arguments.output_folder + "/" + outname
        cmd = "sampler.py -i " + pair[0] + " " + pair[1] + " -c " + str(arguments.coverage_depth)  \
            + " -s " + str(arguments.genome_size) + " -o " + outname
        if arguments.gzip_output:
            cmd += " -g"
        os.system(cmd)
    for single in singles:
        outname = single.split("/")[-1]
        outname = outname.split("_")[0] + "_subsampled" + str(arguments.coverage_depth)
        outname = arguments.output_folder + "/" + outname
        cmd = "sampler.py -i " + single + " -c " + str(arguments.coverage_depth) \
              + " -s " + str(arguments.genome_size) + " -o " + outname
        if arguments.gzip_output:
            cmd += " -g"
        os.system(cmd)



