# Code to subsample fastq files to a desired coverage level.


class Sampler(object):

    def read_fastq(self):
        """
        Reads the fastq files that are input, checking if they're gzipped along the way.
        :return: records - a list of fastq objects.
        """
        from Bio import SeqIO
        records = list()
        try:
            for record in SeqIO.parse(self.fastq_file, "fastq"):
                records.append(record)
        except UnicodeDecodeError:
            import gzip
            fastq = gzip.open(self.fastq_file, 'rt')
            for record in SeqIO.parse(fastq, "fastq"):
                records.append(record)
        return records

    def sample(self):
        """
        Reads in fastq file, then samples reads until desired coverage depth is met. Outputs the reads to self.outfile
        """
        from Bio import SeqIO
        import random
        basepairs_written = 0
        output_list = list()
        print("Reading in fastq file...")
        record_list = self.read_fastq()

        print("Randomly selecting reads...")
        while basepairs_written < self.required_basepairs:
            num = random.randint(0, len(record_list) - 1)
            basepairs_written += len(record_list[num].seq)
            output_list.append(record_list[num])
            del record_list[num]

        print("Writing reads to file...")
        with open(self.outfile, "w") as handle:
            SeqIO.write(output_list, handle, "fastq")

    def gzip_output(self):
        """
        Gzips output, called  when
        :return:
        """
        import gzip
        import os
        print("Gzipping output...")
        f_in = open(self.outfile)
        gzip_out = gzip.open(self.outfile + ".gz", "wt")
        gzip_out.writelines(f_in)
        f_in.close()
        gzip_out.close()
        os.remove(self.outfile)

    def __init__(self, args):
        self.fastq_file = args.input_file
        self.coverage_depth = args.coverage_depth
        self.genome_size = args.genome_size
        self.required_basepairs = self.coverage_depth * self.genome_size
        self.outfile = args.output_file
        self.sample()
        if args.gzip_output:
            self.gzip_output()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="File to be subsampled. Must be in FASTQ format.", required=True)
    parser.add_argument("-c", "--coverage_depth", type=float, required=True, help="Coverage depth to subsample to.")
    parser.add_argument("-s", "--genome_size", type=int, required=True, help="Expected genome size in bp.")
    parser.add_argument("-o", "--output_file", type=str, default="out.fastq", help="Name of your output fastq file."
                                                                                   "Default is out.fastq")
    parser.add_argument("-g", "--gzip_output", default=False, action="store_true", help="When included"
                                                                                        ", gzips output.")
    arguments = parser.parse_args()
    Sampler(arguments)
