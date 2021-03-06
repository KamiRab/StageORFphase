import argparse
import errno
import os
import subprocess
import time
from datetime import datetime
import pandas as pd
import Mapper
from Bam2Reads_function.BAM2Reads import BAM2Reads

R2 = "CTGTAGGCACCATCAAT"  # TODO regarder nombre de reads
R3 = "TGGAATTCTCGGGTGCCAAGG"


# todo maybe just isolate the B2R/plots
try:
    from subprocess import CompletedProcess
except ImportError:
    # Python 2

    class CompletedProcess:

        def __init__(self, args, returncode, stdout=None, stderr=None):
            self.args = args
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

        def check_returncode(self):
            if self.returncode != 0:
                err = subprocess.CalledProcessError(self.returncode, self.args, output=self.stdout)
                raise err
            return self.returncode


    def sp_run(*popenargs, **kwargs):
        input = kwargs.pop("input", None)
        check = kwargs.pop("handle", False)
        if input is not None:
            if 'stdin' in kwargs:
                raise ValueError('stdin and input arguments may not both be used.')
            kwargs['stdin'] = subprocess.PIPE
        process = subprocess.Popen(*popenargs, **kwargs)
        try:
            outs, errs = process.communicate(input)
        except:
            process.kill()
            process.wait()
            raise
        returncode = process.poll()
        if check and returncode:
            raise subprocess.CalledProcessError(returncode, popenargs, output=outs)
        return CompletedProcess(popenargs, returncode, stdout=outs, stderr=errs)


    subprocess.run = sp_run
    # ^ This monkey patch allows it work on Python 2 or 3 the same way

def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORF genetic elements mapper')
    parser.add_argument("-gff",
                        type=str,
                        required=True,
                        nargs="*",
                        help="GFF annotation file")
    parser.add_argument("-fasta",
                        type=str,
                        required=True,
                        nargs="*",
                        help="FASTA file containing the genome sequence")
    parser.add_argument("-fastq",
                        type=str,
                        required=True,
                        nargs="*",
                        help="FASTQ file containing the riboseq reads")
    parser.add_argument("-kmer",
                        type=int,
                        required=True,
                        nargs="*",
                        help="Kmer to analysis")
    parser.add_argument("-cutdir",
                        type=str,
                        required=False,  # if directory not given, launch cutadapt
                        nargs=1,
                        help="Directory containing the directories of reads filtered by size. "
                             "The subdirectories need to be named as kmer_n with n the size of the reads"
                             "If you already launched Ribomap one time with the use of cutadapt activated, "
                             "the directory output should be named like the -o argument"
                        )
    parser.add_argument("-thr",  # can only choose 1 threshold for all the files
                        type=int,
                        required=False,
                        nargs=1,
                        default=10,
                        help="Minimal number of reads to consider to calculate the mean and median per phase"
                        )
    parser.add_argument("-adapt",
                        type=str,
                        required=False,  # if not given, will look in list and if not found error message (TODO)
                        nargs="*",
                        help="Nucleotidic sequence of the adaptor for riboseq")
    parser.add_argument("-features_include",
                        type=str,
                        required=False,
                        nargs="*",
                        default=["all"],
                        help="Type(s) of genetic elements to consider. By default it is all "
                             "For all you can enter \"all\" and for just non coding elements \"nc\""
                             "Warning : if you enter multiple files, all the elements will be included for each file")
    parser.add_argument("-features_exclude",
                        type=str,
                        required=False,
                        nargs="*",
                        default=["None"],
                        help="Type(s) of genetic elements to not consider. By default it is none "
                             "Warning : if you enter multiple files, all the elements will be exluded for each file")
    parser.add_argument("-options",
                        type=str,
                        nargs="?",
                        default="IMBP",
                        help='''What needs to be launched:
                            I : indexing of the genome file (Bowtie-build)
                            M : mapping of the reads on the index file to bam files
                            B : counting of the reads by phase and periodicity (Bam2reads function)
                            P : plotting of periodicity and phasing
                            WARNING : I,M and B should be keeped if you never launched Ribomap or the name/path 
                            of the files, may differ from what is generated by Ribomap. 
                        ''')
    parser.add_argument("-o",
                        type=str,
                        # action='store',
                        required=True,
                        nargs="?",
                        help="Output name of analysis file to put after the genome name. "
                             "For example, you can put \"CDS\" to have Ecoli_CDS. If you didn't enter a directory for "
                             "the cutadapt files, the directory name of the cutadapt files will be this")
    args = parser.parse_args()
    return args


def main():
    start_time = datetime.now()

    # Get the parameters
    parameters = get_args()
    Mapper.parameters_verification(parameters)
    typein = parameters.features_include
    typeout = parameters.features_exclude
    kmer = parameters.kmer
    output_name = parameters.o
    reads_thr = parameters.thr
    options = list(parameters.options)
    for input_file in range(len(parameters.fastq)):
        start_time_input =datetime.now()
        if len(parameters.fasta) == 1:
            genome_file = parameters.fasta[0] #path to genome file
        else:
            genome_file = parameters.fasta[input_file]
        genome_name = os.path.basename(genome_file)
        genome_name = os.path.splitext(genome_name)[0] #name of the genome

        riboseq_file = parameters.fastq[input_file] #path of the fastq file
        riboseq_name = os.path.basename(riboseq_file)
        riboseq_name = os.path.splitext(riboseq_name)[0] #name of the riboseq experience

        if len(parameters.gff) == 1:
            gff_file = parameters.gff[0]
        else:
            gff_file = parameters.gff[input_file]
        with open("ORFribomap_{}_{}_log.txt".format(riboseq_name,output_name), "w") as log:
            # 1. Annotation of the ORFs with ORFtrack
            mapping_gff = "mapping_orf_{}.gff".format(genome_name)
            if os.path.exists(mapping_gff): #if the ORFtrack file exists already, ORFtrack is skipped
                print(mapping_gff + " exists, it will be used to extract the intergenic ORFs. If the file doesn't result "
                                    "from orftrack, please rename it and launch ORFribomap again".format(genome_name))
            else:
                orftrack_cmd = "orftrack -fna {} -gff {}".format(genome_file, gff_file)
                print("Launch : ", orftrack_cmd)
                process_orftrack = subprocess.run(orftrack_cmd, shell=True, universal_newlines=True)

            # 3. Potential launch of cutadapt and definition of the directory that will contain all the analysis files
            if not parameters.cutdir:
                cutdir = "."
                if len(parameters.adapt) == 1: #the adaptor is used for all the files
                    adapt = parameters.adapt[0]
                else:
                    adapt = parameters.adapt[input_file] #each file has its adaptor
                print(" The cutadapt process will be launched")
                try:
                    print("\nNo directory for cutadapt files has been given, the cutadapt process is launching. "
                          "The reads are cut in {}-kmers".format(parameters.kmer))
                    for size in kmer:
                        log.write(Mapper.cut_reads(size, riboseq_file, adapt, riboseq_name))
                except:
                    print('''cutadapt seems to not be installed''')
            else: #there exist a directory with the the reads already trimmed and cut
                print("\nYou entered directory(ies) for cutadapt files. No cutadapt process will be launched")
                time.sleep(1)
                if len(parameters.cutdir) == 1:
                    cutdir = parameters.cutdir[0]
                else :
                    cutdir = parameters.cutdir[input_file]
                for size in kmer:
                    name_dir = "{}/kmer_{}".format(cutdir, size)
                    if not os.path.isdir(name_dir):
                        print(name_dir + " does not exist. The directories need to be named as kmer_x with x the size of "
                                         "the kmer. ")
                        exit()

            # 4.Mapping of the reads on the non-translated sequences
            print("We start the mapping")

            # a. Building of the Bowtie index from the orfget output
            if "I" in options:
                cmd_bowtie = 'bowtie-build {} {}_all'.format(genome_file, genome_name)
                print("Command launched: ", cmd_bowtie)
                process_bowtie = subprocess.run(cmd_bowtie, shell=True, universal_newlines=True, stdout=subprocess.PIPE)
                log.write(process_bowtie.stdout)
            else:
                print("The basename of Bowtie index files are expected to be : {}_all".format(genome_name))

            # b->f. Mapping of the reads on the Bowtie index
            if "M" in options:
                try: #verification that pysam and bokeh are installed
                    import pysam
                    import bokeh
                except ImportError:
                    print("You need to install the python packages pysam and bokeh")
                for size in kmer: #mapping of the reads
                    log.write(Mapper.map2bam(cutdir,size, genome_name, riboseq_name, "all"))
            else:
                print("The name and path of the bam files are expected to be :"
                      "{}/kmer_n/{}_kmer_n_all_sorted_mapped.bam with n the size of the reads".format(cutdir, riboseq_name))

            # g. Creation of the counting tables and periodicity tables
            if "B" in options:
                print("Counting tables are generated")
                log.write(BAM2Reads(riboseq_name, mapping_gff, kmer, output_name, typein, typeout))
            else:
                print("The name and path of the counting and periodicity files are expected to be :\n"
                      "{}/kmer_n/{}_kmer_n_{}_reads.tab with n the size of the reads\n"
                      "{}/kmer_n/{}_kmer_n_{}_periodicity_start.tab with n the size of the reads\n"
                      "{}/kmer_n/{}_kmer_n_{}_periodicity_stop.tab with n the size of the reads\n".format(cutdir,
                                                                                                          riboseq_name,
                                                                                                          output_name,
                                                                                                          cutdir,
                                                                                                          riboseq_name,
                                                                                                          output_name,
                                                                                                          cutdir,
                                                                                                          riboseq_name,
                                                                                                          output_name))

            # 5. Plotting of phasing, periodicity of the start and the stop
            if "P" in options:
                for size in kmer:
                    tab = pd.read_table(
                        "{}/kmer_{}/{}_kmer_{}_{}_reads.tab".format(cutdir, size, riboseq_name, size, output_name),
                        sep='\t')
                    # Plot of the reads phase and periodicity
                    Mapper.reads_phase_plot(tab, size, riboseq_name,
                                            reads_thr, output_name)
                    Mapper.reads_periodicity(size, riboseq_name, output_name, "start")
                    Mapper.reads_periodicity(size, riboseq_name, output_name, "stop")
            end_time_input = datetime.now()
            log.write("Duration for {} :{}".format(riboseq_name, end_time_input-start_time_input))
    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))


main()
