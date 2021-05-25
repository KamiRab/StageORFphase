#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from datetime import datetime


def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORF phasing of genetic elements')
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
    parser.add_argument("-kmer",  # is valable for all files
                        type=int,
                        required=False,
                        nargs=2,
                        default=[26, 30],
                        help="Range of size of reads that will be tested to find the size with maximal phase 0")
    parser.add_argument("-cutdir",
                        type=str,
                        required=False,  # if directory not given, launch cutadapt
                        nargs=1,
                        help="Directory containing the directories of cutadapt files"
                        )
    parser.add_argument("-adapt",
                        type=str,
                        required=False,  # if not given, will look in list and if not found error message (TODO)
                        nargs="*",
                        help="Nucleotidic sequence of the adaptor for riboseq")
    parser.add_argument("-typein",
                        type=str,
                        required=False,
                        nargs="*",
                        default="CDS",
                        help="Type of genetic elements to consider with ORFmine nomenclature. "
                             "For all you can enter \"all\" and for just non coding elements \"nc\""
                             "Warning : if you enter multiple files, all the elements will be included for each file")
    parser.add_argument("-typeout",
                        type=str,
                        required=False,
                        nargs="*",
                        default="CDS",
                        help="Type of genetic elements to not consider with ORFmine nomenclature. "
                             "Warning : if you enter multiple files, all the elements will be exluded for each file")
    parser.add_argument("-")
    parser.add_argument("-o",
                        type=str,
                        #action='store',
                        required=True,
                        nargs=1,
                        help="Output name of analysis file to put after the genome name. "
                             "For example, you can put \"CDS\" to have Ecoli_CDS")
    parser.add_argument("-custom", #TODO options to declare
                        type=list,
                        required=False,
                        nargs="?",
                        default=["PR"],
                        help=
                        '''What to test :
                            P : search best reads size
                            R : map ORFs found by ORFtrack on genome
                        ''')
    args = parser.parse_args()
    return args

def main():
    start_time = datetime.now()

    # Get the parameters
    global parameters
    parameters = get_args()
    Mapper.parameters_verification(parameters)
    typein = parameters.typein
    typeout = parameters.typeout
    kmer = parameters.kmer
    output_name = parameters.o
    for input_file in range(len(parameters.gff)):
        genome_file = parameters.fasta[input_file]
        genome_name = os.path.basename(genome_file)
        genome_name = os.path.splitext(genome_name)[0]

        riboseq_file = parameters.fastq[input_file]
        riboseq_name = os.path.basename(riboseq_file)
        riboseq_name = os.path.splitext(riboseq_name)[0]

        gff_file = parameters.gff[input_file]

        # 1. Annotation of the ORFs with ORFtrack
        mapping_gff = "mapping_orf_{}.gff".format(genome_name)
        if os.path.exists(mapping_gff):
            print(mapping_gff + " exists, it will be used to extract the intergenic ORFs. "
                  "If the file doesn't result from orftrack, please rename it and launch ORFphase again".format(
                genome_name))
        else:
            orftrack_cmd = "orftrack -fna {} -gff {}".format(genome_file, gff_file)
            print("Launch : ", orftrack_cmd)
            process_orftrack = subprocess.run(orftrack_cmd, shell=True)

        # 2. Extraction of the non-translated sequences of ORF:
        elements_in = " ".join(typein)
        elements_out = " ".join(typeout)

        orfget_cmd = "orfget -fna {} -gff {} -features_include {} -features_exclude {} -o {}_{} -type nucl".format(
            genome_file, mapping_gff, elements_in, elements_out, genome_name, output_name)
        print("Launch :\t", orfget_cmd)
        process_orfget = subprocess.run(orfget_cmd, shell=True)

        # 3. Potential launch of cutadapt and definition of cutadapt files directory
        if not parameters.cutdir:
            if len(parameters.adapt) == 1:
                adapt = parameters.adapt[0]
            adapt = parameters.adapt[input_file]
            print(" The cutadapt process will be launched")
            try:
                import cutadapt
                print("\nNo directory for cutadapt files has been given, the cutadapt process is launching. "
                      "The reads are cut in {}-kmers".format(parameters.kmer))
            except ImportError:
                print('''cutadapt is not installed, please check ''')  # add link for installation
            cutdir = "./" + output_name
            try:
                os.mkdir(cutdir)  # Create nc_intergenic directory to distinguish them from potential CDS analysis
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
            for i in kmer:
                Mapper.cut_reads(i, riboseq_file, adapt, riboseq_name, cutdir)
        else:
            print("\nYou entered directory(ies) for cutadapt files. No cutadapt process will be launched")
            cutdir = parameters.cutdir[input_file]
            for i in kmer:
                name_dir = "{}/kmer_{}".format(cutdir, i)
                if not os.path.isdir(name_dir):
                    print("The directories need to be named as kmer_x with x the size of the kmer.")
                    exit()

        # 4.Mapping of the reads on the non-translated sequences
        print("We start the mapping")

        # a. Building of the Bowtie index from the orfget output
        cmd_bowtie = 'bowtie-build {}_{}.nfasta {}_{}'.format(genome_name, output_name, genome_name, output_name)
        process_bowtie = subprocess.run(cmd_bowtie, shell=True)

        # b->g. Create count table
        try:
            import pysam
            import bokeh
        except ImportError:
            print("You need to install the python packages pysam and bokeh")
        for i in kmer:
            Mapper.map2bam(i, genome_name, riboseq_name, cutdir, parameters.type, gff_file)

        # 5. Find best kmer cut to have phase 0
        p0_by_kmer = {}
        best_kmer = {}
        for i in kmer:
            tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_reads.tab".format(cutdir, kmer, riboseq_name, kmer), sep='\t',
                                names=["ID", "Number of reads", "Number of p0", "Number of p1", "Number of p2",
                                       "Percentage of p0", "Percentage of p1", "Percentage of p2"])
            # Plot of the reads phase and periodicity
            Mapper.reads_phase_plot(tab, kmer, riboseq_name,
                                    10)  # for CDS, we take in consideration only genes with more than 100 reads
            Mapper.reads_periodicity(kmer, riboseq_name, cutdir, "start")
            Mapper.reads_periodicity(kmer, riboseq_name, cutdir, "stop")

    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))

