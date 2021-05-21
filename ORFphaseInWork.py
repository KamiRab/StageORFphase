#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import concurrent.futures.process
import errno
import os
import subprocess
import pandas as pd
from datetime import datetime
import Mapper
import detect_igr.BAM2Reads


# TODO R2 with good adapter
#  TODO ORFget + extend (parametre) ==> CDS ==> periodicite
# detectIGR : pysam, bokeh
# TODO partial codon ?
# install cutadapt by command "pip3 install cutadapt"
# install bowtie by https://github.com/BenLangmead/bowtie then make then "sudo make install"
# install samtools through downloading the current source release of samtools (www.htslib.org/download/)
def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORFphase')
    parser.add_argument("-gff",
                        type=str,
                        required=True,
                        nargs="+",
                        help="GFF annotation file")
    parser.add_argument("-fasta",
                        type=str,
                        required=True,
                        nargs="+",
                        help="FASTA file containing the genome sequence")
    parser.add_argument("-fastq",
                        type=str,
                        required=True,
                        nargs="+",
                        help="FASTQ file containing the riboseq reads")
    parser.add_argument("-kmer",  # can only choose one option of kmer
                        type=int,
                        required=False,
                        nargs=2,
                        default=[26, 30])
    parser.add_argument("-cutdir",
                        type=str,
                        required=False,  # if directory not given, launch cutadapt
                        nargs="+",
                        help="Directory containing the directories of cutadapt files"
                        )
    parser.add_argument("-adapt",
                        type=str,
                        required=False,  # if not given, will look in list and if not found error message (TODO)
                        nargs="+",
                        help="Nucleotidic sequence of the adaptor for riboseq")
    parser.add_argument("-thr",  # can only choose 1 threshold for all the files
                        type=int,
                        required=False,
                        choices=range(0, 100),
                        nargs=1,
                        default=90,
                        help="Threshold for keeping the reads in phase 0")
    parser.add_argument("-type",
                        type=str,
                        required=False,
                        choices=["CDS", "nc_intergenic"],
                        nargs="?",
                        default="CDS",
                        help="Analysis on CDS (CDS) or intergenic reads (nc_intergenic)")
    args = parser.parse_args()
    return args


# 2. Create gff file of the transcriptome or intergenic ORFs
def read_multiFASTA(fasta_file):
    """
    Create dictionary from fasta_file with the gene name as key and the value being the sequence
    :param fasta_file: fasta file to read
    :return: dictionnary
    """
    dico = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                name = str(line.split()[0])[1:]
                dico[name] = ''
            elif line == '\n':
                continue
            else:
                seq = line.strip().replace("*", "")
                dico[name] = dico[name] + seq
    return (dico)


def reads_phase_percentage(tab):
    """
    From the reads count table, give the percentage of all phases
    :param tab: reads count table
    :return: list of the percentage of the three phases
    """
    perc = []
    # Sum of phase columns in tab
    phase0sum = tab["Number of p0"].sum()
    phase1sum = tab["Number of p1"].sum()
    phase2sum = tab["Number of p2"].sum()
    sumreads = tab["Number of reads"].sum()

    # List with percentage of each phase
    perc.append(phase0sum / sumreads)
    perc.append(phase1sum / sumreads)
    perc.append(phase2sum / sumreads)
    return perc


def main():
    start_time = datetime.now()

    # Get the parameters
    global parameters
    parameters = get_args()
    Mapper.parameters_verification(parameters)
    for input_file in range(len(parameters.gff)):
        kmer = range(parameters.kmer[0], parameters.kmer[1] + 1)

        genome_file = parameters.fasta[input_file]
        genome_name = os.path.basename(genome_file)
        genome_name = os.path.splitext(genome_name)[0]

        riboseq_file = parameters.fastq[input_file]
        riboseq_name = os.path.basename(riboseq_file)
        riboseq_name = os.path.splitext(riboseq_name)[0]

        gff_file = parameters.gff[input_file]

        # 1. Extraction of the non-translated sequences of CDS:
        orfget_cmd = "orfget -fna {} -gff {} -features_include CDS -o {}_CDS -type nucl".format(
            genome_file, gff_file, genome_name)
        print("Extract the Transcriptome:")
        print("Launch :\t", orfget_cmd)
        process_orfget = subprocess.run(orfget_cmd, shell=True)

        # 2. Generation of pseudo GFF file of the CDS transcriptome
        gff_to_compare = genome_name + "_transcriptome.gff"
        transcriptome = read_multiFASTA(fasta_file=genome_name + "_CDS.nfasta")
        with open(gff_to_compare, "w") as wgff:
            for i in transcriptome:
                wgff.write('{:20s}\t{}\t{:20s}\t{:d}\t{:d}\t{}\t{}\t{}\t{}\n'.format(
                    i, 'SGD', 'gene', 1, len(transcriptome[i]), '.', '+', '0', str('ID=' + i)))

        # 3. Potential launch of cutadapt and definition of cutadapt files directory
        if not parameters.cutdir:
            cutdir = "."
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
            for i in kmer:
                Mapper.cut_reads(i, riboseq_file, adapt, riboseq_name, cutdir)
        else:
            print("\nYou entered directory(ies) for cutadapt files. No cutadapt process will be launched")
            if len(parameters.cutdir) == 1:
                cutdir = parameters.cutdir[0]
            else :
                cutdir = parameters.cutdir[input_file]
            for i in kmer:
                name_dir = "{}/kmer_{}".format(cutdir, i)
                if not os.path.isdir(name_dir):
                    print("The directories need to be named as kmer_x with x the size of the kmer.")
                    exit()

        # 4.Mapping of the reads on the non-translated sequences
        print("We start the mapping")

        # a. Building of the index
        cmd_bowtie = 'bowtie-build {}_CDS.nfasta {}_CDS'.format(genome_name, genome_name)
        process_bowtie = subprocess.run(cmd_bowtie, shell=True)

        # b->g. Create count table
        try:
            import pysam
            import bokeh
        except ImportError:
            print("You need to install the python packages pysam and bokeh")
        for i in kmer:
            Mapper.map_in_bam_and_count(i, genome_name, riboseq_name, cutdir, "CDS", gff_to_compare)
        # with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
        #     executor.map(map_in_bam_and_count, kmer, [gname] * len(kmer), [rname] * len(kmer), [cutdir] * len(kmer),
        #     [parameters.type]*len(kmer),[gff_to_compare]*len(kmer))

        # 5. Find best kmer cut to have phase 0
        p0_by_kmer = {}
        best_kmer = {}
        for i in kmer:
            tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_reads.tab".format(cutdir, i, riboseq_name, i), sep='\t',
                                names=["ID", "Number of reads", "Number of p0", "Number of p1", "Number of p2",
                                       "Percentage of p0", "Percentage of p1", "Percentage of p2"])
            # Plot of the reads phase and periodicity
            Mapper.reads_phase_plot(tab, i, riboseq_name,
                                    100)  # for CDS, we take in consideration only genes with more than 100 reads
            Mapper.reads_periodicity(i, riboseq_name, cutdir, "start")
            Mapper.reads_periodicity(i, riboseq_name, cutdir, "stop")
            p0_by_kmer[i] = reads_phase_percentage(tab)
            print(
                "For the {}_kmer_{}, there are : \n{:.3f}% reads in phase 0 \n{:.3f}% reads in phase 1 \n"
                "{:.3f}% reads in phase 2".format(riboseq_name, i, p0_by_kmer[i][0],
                                                  p0_by_kmer[i][1], p0_by_kmer[i][2]))
            if p0_by_kmer[i][0] > parameters.thr / 100:
                best_kmer[i] = p0_by_kmer[i][0]
        print("For {} the kmers with a mean superior to {}% are: ".format(riboseq_file, parameters.thr), best_kmer)
        if not best_kmer:
            print("No kmer had a mean of phase 0 superior to the threshold for {}".format(riboseq_file))
        # kmer_choosen = input(
        #         "Kmer(s) to use (if you input multiple kmer, put a space between them)")
        #     list_kmer_chooser = kmer_choosen.split()
        #     # suppress_other_directories = input("Do you want to suppress the cutadapt directories for the non chooser "
        #                                    "kmer ? [Y/N]")
        # if suppress_other_directories == "Y" or suppress_other_directories == "y":
        #     for i in kmer:
        #         if i not in kmer_choosen:
        #             os.rmdir("kmer_{}".format(i))

    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))


main()
