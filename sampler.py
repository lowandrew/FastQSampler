#!/usr/bin/env python


class Sampler(object):

    @staticmethod
    def read_fastq(fastq_file):
        """
        Reads the fastq files that are input, checking if they're gzipped along the way.
        :return: records - a list of fastq objects.
        """
        from Bio import SeqIO
        # Empty list that will hold my fastqs
        records = list()
        # Use SeqIO to read in all my records.
        try:
            for record in SeqIO.parse(fastq_file, "fastq"):
                records.append(record)
        # For some reason different gzipped files give different exceptions, so we'll just except everything.
        except:
            import gzip
            fastq = gzip.open(fastq_file, 'rt')
            for record in SeqIO.parse(fastq, "fastq"):
                records.append(record)
        return records

    def sample_single_reads(self):
        """
        Reads in fastq file, then samples reads until desired coverage depth is met. Outputs the reads to self.outfile
        """
        from Bio import SeqIO
        import random
        import sys
        # Keep track of how many basepairs we've written so we know when to stop sampling.
        basepairs_written = 0
        output_list = list()
        # Use our wonderful read_fastq method to get a list of all our fastq seqs.
        print("Reading in fastq file " + self.fastq_file[0] + "...")
        record_list = Sampler.read_fastq(self.fastq_file[0])

        print("Randomly selecting reads...")
        # Sample reads until we've gotten enough coverage.
        while basepairs_written < self.required_basepairs:
            # It's possible people will try to sample too deeply (i.e. sample to 50X when they only have 30X coverage.
            # If they do, this try/except will boot them with an appropriate error message.
            try:
                num = random.randint(0, len(record_list) - 1)
            except ValueError:
                print("It looks like you're trying to sample reads to a depth that is too deep! Exiting...")
                sys.exit()
            # Add to the running total of base pairs sampled.
            basepairs_written += len(record_list[num].seq)
            output_list.append(record_list[num])
            # Delete the record we just sampled so it doesn't get re-sampled.
            del record_list[num]

        # Write reads to file.
        print("Writing reads to file...")
        with open(self.outfile + ".fastq", "w") as handle:
            SeqIO.write(output_list, handle, "fastq")

        if self.gzip_output:
            Sampler.gzip_output(self.outfile + ".fastq")

    def sample_paired_reads(self):
        """
        Reads in fastq file, then samples reads until desired coverage depth is met. Outputs the reads to self.outfile
        """
        from Bio import SeqIO
        import random
        import sys
        # Keep track of how many basepairs we've written so we know when to stop sampling.
        basepairs_written = 0
        output_list_forward = list()
        output_list_reverse = list()
        # Use our wonderful read_fastq method to get a list of all our fastq seqs.
        print("Reading in fastq file " + self.fastq_file[0] + "...")
        record_list_forward = Sampler.read_fastq(self.fastq_file[0])
        print("Reading in fastq file " + self.fastq_file[1] + "...")
        record_list_reverse = Sampler.read_fastq(self.fastq_file[1])

        print("Randomly selecting reads...")
        # Sample reads until we've gotten enough coverage.
        while basepairs_written < self.required_basepairs:
            # It's possible people will try to sample too deeply (i.e. sample to 50X when they only have 30X coverage.
            # If they do, this try/except will boot them with an appropriate error message.
            try:
                num = random.randint(0, len(record_list_forward) - 1)
            except ValueError:
                print("It looks like you're trying to sample reads to a depth that is too deep! Exiting...")
                sys.exit()
            # Add to the running total of base pairs sampled. Multiply by wo since we're in paired end mode.
            basepairs_written += len(record_list_forward[num].seq) * 2
            output_list_forward.append(record_list_forward[num])
            output_list_reverse.append(record_list_reverse[num])
            # Delete the record we just sampled so it doesn't get re-sampled.
            del record_list_forward[num]
            del record_list_reverse[num]

        # Write reads to file.
        print("Writing forward reads to file...")
        with open(self.outfile + "_R1.fastq", "w") as handle:
            SeqIO.write(output_list_forward, handle, "fastq")
        print("Writing reverse reads to file...")
        with open(self.outfile + "_R2.fastq", "w") as handle:
            SeqIO.write(output_list_reverse, handle, "fastq")

        if self.gzip_output:
            Sampler.gzip_output(self.outfile + "_R1.fastq")
            Sampler.gzip_output(self.outfile + "_R2.fastq")

    @staticmethod
    def gzip_output(output_fastq_file):
        """
        Gzips output, called when user adds the -g flag when calling the program.
        """
        import gzip
        import os
        print("Gzipping " + output_fastq_file + "...")
        f_in = open(output_fastq_file)
        gzip_out = gzip.open(output_fastq_file + ".gz", "wt")
        gzip_out.writelines(f_in)
        f_in.close()
        gzip_out.close()
        os.remove(output_fastq_file)

    def __init__(self, args):
        import sys
        self.fastq_file = args.input_file
        self.coverage_depth = args.coverage_depth
        self.genome_size = args.genome_size
        self.required_basepairs = self.coverage_depth * self.genome_size
        self.outfile = args.output_file
        self.gzip_output = args.gzip_output
        if len(self.fastq_file) == 1:
            self.sample_single_reads()
        elif len(self.fastq_file) == 2:
            self.sample_paired_reads()
        else:
            print("Input reads don't appear to be single or paired! Exiting...")
            sys.exit()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Reads to be subsampled. Must be in FASTQ format.",
                        nargs="*", required=True)
    parser.add_argument("-c", "--coverage_depth", type=float, required=True, help="Coverage depth to subsample to.")
    parser.add_argument("-s", "--genome_size", type=int, required=True, help="Expected genome size in bp.")
    parser.add_argument("-o", "--output_file", type=str, default="out", help="Base name of your output fastq file. "
                                                                             "Default is out")
    parser.add_argument("-g", "--gzip_output", default=False, action="store_true", help="When included"
                                                                                        ", gzips output.")
    arguments = parser.parse_args()
    Sampler(arguments)
