import argparse
import concurrent.futures.process
import errno
import os
import subprocess
import time
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


# TODO R2 strange
# detectIGR : pysam, bokeh
# TODO partial codon ?
def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORF phase')
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
    parser.add_argument("-kmer",  # can only choose one option of kmer
                        type=int,
                        required=False,
                        nargs=2,
                        default=[26, 30])
    parser.add_argument("-cutdir",
                        type=str,
                        required=False,  # if directory not given, launch cutadapt
                        nargs="*",
                        help="Directory containing the directories of cutadapt files"
                        )
    parser.add_argument("-adapt",
                        type=str,
                        required=False,  # if not given, will look in list and if not found error message (TODO)
                        nargs="*",
                        help="Nucleotidic sequence of the adaptor for riboseq")
    parser.add_argument("-thr",  # can only choose 1 threshold for all the files
                        type=int,
                        required=False,
                        nargs=1,
                        default=90,
                        help="Threshold for keeping the reads in phase 0")
    args = parser.parse_args()
    return args


def file_verification(parameters):
    gff = parameters.gff
    fasta = parameters.fasta
    adapt = parameters.adapt
    fastq = parameters.fastq
    if len(gff) == len(fasta) and len(fasta) == len(fastq) and len(fastq) == len(adapt):
        print("The associations of the files are:")
        for i in range(len(gff)):
            print("GFF: {} \tFASTA: {}\tFASTAQ: {}\t adaptor: {}".format(gff[i], fasta[i], fastq[i], adapt[i]))
    else:
        print("There is something wrong with the number of files input. There must be the same number of gff, "
              "fasta files, adaptors sequences and fastq files")
        exit()


def classify_reads(kmer, fastq_file, adaptor, rname, cutdir):
    '''

    :param kmer:
    :param fastq_file:
    :param adaptor:
    :param rname:
    :param cutdir:
    :return:
    '''
    try:
        os.mkdir("kmer_{}".format(str(kmer)))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    cutadapt_cmd = "cutadapt {} -j 1 --quality-base=33 -a {} \
-m {} -M {} -e 0.12 -o {}/kmer_{}/{}_kmer{}.fastq".format(fastq_file, adaptor, str(kmer), str(kmer), cutdir,
                                                          str(kmer),
                                                          rname, str(kmer))
    print("We are cutting the reads in {}-mers".format(kmer))
    process = subprocess.run(cutadapt_cmd, shell=True)


# 2. Create gff file of the transcriptome
def read_multiFASTA(fasta_file):
    '''

    :param fasta_file:
    :return:
    '''
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


def map_in_bam_and_count(kmer, gname, rname, cutdir):
    '''

    :param kmer:
    :param gname:
    :param rname:
    :param cutdir:
    :return:
    '''

    # Commands
    # b. Map the reads on the transcriptome and generate a sam file
    cmd_bowtie = "bowtie --wrapper basic-0 -v 2 -y -a -m 1 --best --strata -p 18 -x {}_CDS {}/kmer_{}/{}_kmer{}.fastq -S kmer_{}/{}_kmer{}.sam".format(
        gname, cutdir, str(kmer), rname, str(kmer), str(kmer), rname, str(kmer))
    # c. Transform sam into bam
    cmd_sam_bam = "samtools view -h -b -S {}/kmer_{}/{}_kmer{}.sam > {}/kmer_{}/{}_kmer{}.bam".format(cutdir, str(kmer),
                                                                                                      rname,
                                                                                                      str(kmer), cutdir,
                                                                                                      str(kmer), rname,
                                                                                                      str(kmer))
    # d. Sort the bam file
    cmd_sam_sort = "samtools sort {}/kmer_{}/{}_kmer{}.bam -o {}/kmer_{}/{}_kmer{}_sorted.bam".format(cutdir, str(kmer),
                                                                                                      rname,
                                                                                                      str(kmer), cutdir,
                                                                                                      str(kmer), rname,
                                                                                                      str(kmer))
    # e. Keep only the mapped reads
    cmd_sam_map = "samtools view -F 4 -f 0,16 -h -b {}/kmer_{}/{}_kmer{}_sorted.bam > {}/kmer_{}/{}_kmer{}_sorted_mapped.bam".format(
        cutdir, str(kmer), rname, str(kmer), cutdir, str(kmer), rname, str(kmer))
    # f. Index the bam file
    cmd_sam_index = "samtools index {}/kmer_{}/{}_kmer{}_sorted_mapped.bam".format(cutdir, str(kmer), rname, str(kmer))

    # g. generate count table
    cmd_count = "python ./detect_igr/BAM2Reads.py -bam kmer_{}/{}_kmer{}_sorted_mapped.bam -gff {}_transcriptome.gff " \
                "-output_path kmer_{} -output_name {}_kmer_{}".format(
        str(kmer), str(rname), str(kmer), str(gname), str(kmer), str(rname), str(kmer))

    # Processes
    print("Command launched: ", cmd_bowtie)
    process1 = subprocess.run(cmd_bowtie, shell=True)

    print("Command launched: ", cmd_sam_bam)
    process2 = subprocess.run(cmd_sam_bam, shell=True)

    print("Command launched: ", cmd_sam_sort)
    process3 = subprocess.run(cmd_sam_sort, shell=True)

    print("Command launched: ", cmd_sam_map)
    process4 = subprocess.run(cmd_sam_map, shell=True)

    print("Command launched: ", cmd_sam_index)
    process5 = subprocess.run(cmd_sam_index, shell=True)

    print("Command launched: ", cmd_count)
    process6 = subprocess.run(cmd_count, shell=True)

    print("The count table has been generated for kmer {}".format(kmer))


def mean_phase(kmer, rname, cutdir):
    '''

    :param kmer:
    :param rname:
    :param cutdir:
    :return:
    '''
    tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_reads.tab".format(cutdir, kmer, rname, kmer), sep='\t',
                        names=["ID", "Number of reads", "Number of p0", "Number of p1", "Number of p2",
                               "Percentage of p0", "Percentage of p1", "Percentage of p2"])
    perc = []
    # Sum of phase columns in tab
    phase0sum = tab["Number of p0"].sum()
    phase1sum = tab["Number of p1"].sum()
    phase2sum = tab["Number of p2"].sum()
    sumreads = phase0sum + phase1sum + phase2sum
    sumreads2 = tab["Number of reads"].sum()

    # List with frequency of each phase
    perc.append(phase0sum / sumreads)
    perc.append(phase1sum / sumreads)
    perc.append(phase2sum / sumreads)
    print(
        "For the {}_kmer_{}, there are : \n{:.3f}% reads in phase 0 \n{:.3f}% reads in phase 1 \n{:.3f}% reads in phase 2".format(
            rname, kmer, perc[0], perc[1], perc[2]))
    return perc

def distribution_phase_plot(kmer, rname, cutdir, reads_thr):
    '''

    :param kmer:
    :param rname:
    :param cutdir:
    :param reads_thr:
    :return:
    '''
    tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_reads.tab".format(cutdir, kmer, rname, kmer), sep='\t',
                        names=["ID", "Number of reads", "Number of p0", "Number of p1", "Number of p2",
                               "Percentage of p0", "Percentage of p1", "Percentage of p2"])
    tab_select = tab[tab['Number of reads']>reads_thr]
    plt.figure()
    tab_select.boxplot(column=["Percentage of p0", "Percentage of p1", "Percentage of p2"])
    # plt.show()
    plt.savefig('./Boxplot_phases_{}_kmer{}.png'.format(rname, kmer))


def main():
    start_time = datetime.now()

    # Get the parameters
    global parameters
    parameters = get_args()
    file_verification(parameters)
    for input_file in range(len(parameters.gff)):
        kmer = range(parameters.kmer[0], parameters.kmer[1] + 1)

        genome_file = parameters.fasta[input_file]
        gname = os.path.basename(genome_file)
        gname = os.path.splitext(gname)[0]

        fastq_file = parameters.fastq[input_file]
        rname = os.path.basename(fastq_file)
        rname = os.path.splitext(rname)[0]

        gff_file = parameters.gff[input_file]
        adapt = parameters.adapt[input_file]


        # 1. Extraction of the NT sequences of the CDS
        orfget_cmd = "orfget -fna {} -gff {} -features_include CDS -o {}_CDS -type nucl".format(
            genome_file, gff_file, gname)
        print("Extract the Transcriptome:")
        print("Launch :|t", orfget_cmd)
        process_orfget = subprocess.run(orfget_cmd, shell=True)

        # 2. Generation of pseudo GFF file of the CDS transcriptome
        transcriptome = read_multiFASTA(fasta_file=gname + "_CDS.nfasta")
        with open(gname + "_transcriptome.gff", "w") as wgff:
            for i in transcriptome:
                wgff.write('{:20s}\t{}\t{:20s}\t{:d}\t{:d}\t{}\t{}\t{}\t{}\n'.format(
                    i, 'SGD', 'gene', 1, len(transcriptome[i]), '.', '+', '0', str('ID=' + i)))

        # 3. Potential launch of cutadapt and definition of cutadapt files directory
        if not parameters.cutdir:
            cutdir = "."
            print(" The cutadapt process will be launched")
            try:
                import cutadapt
                print("\nNo directory for cutadapt files has been given, the cutadapt process is launching. "
                      "The reads are cut in {}-kmers".format(parameters.kmer))
            except ImportError:
                print('''cutadapt is not installed, please check ''')  # add link for installation

            # for i in kmer:
            #     classify_reads(i, fastq_file, adapt, rname, cutdir)
            with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
                executor.map(classify_reads, kmer, [fastq_file] * len(kmer), [adapt] * len(kmer),
                             [rname] * len(kmer),[cutdir]*len(kmer))
        else:
            print("You entered directory(ies) for cutadapt files. No cutadapt process will be launched")
            cutdir = parameters.cutdir[input_file]
            for i in kmer:
                name_dir = "{}/kmer_{}".format(cutdir, i)
                if not os.path.isdir(name_dir):
                    print("The directories need to be named as kmer_x with x the size of the kmer.")
                    exit()

        # 4.Mapping of the reads on the CDS NT sequences
        print("We start the mapping")

        # a. Building of the index
        cmd_bowtie = 'bowtie-build {}_CDS.nfasta {}_CDS'.format(gname, gname)
        process_bowtie = subprocess.run(cmd_bowtie, shell=True)
        while not os.path.exists(gname + "_CDS.1.ebwt"):  # Why ?
            time.sleep(1)

        # b->g. Create count table
        try :
            import pysam
            import bokeh
        except ImportError:
            print("You need to install the python packages pysam and bokeh")
        for i in kmer:
            map_in_bam_and_count(i, gname, rname, cutdir)
        # with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
        #     executor.map(map_in_bam_and_count, kmer, [gname] * len(kmer), [rname] * len(kmer), [cutdir] * len(kmer))

        # 5. Find best kmer cut to have phase 0
        p0_by_kmer = {}
        best_kmer = {}
        for i in kmer:
            p0_by_kmer[i] = mean_phase(i, rname, cutdir)
            if p0_by_kmer[i][0] > parameters.thr / 100:
                best_kmer[i] = p0_by_kmer[i][0]
            distribution_phase_plot(kmer, rname, cutdir, 100)
        print("The kmer with a mean superior to {}% are: ".format(parameters.thr), best_kmer)
        if best_kmer:
            kmer_choosen = int(input("Kmer to use (Please enter just the number)"))  # should be list (TODO)
            suppress_other_directories = input("Do you want to suppress the cutadapt directories for the non chooser "
                                               "kmer ? [Y/N]")
            if suppress_other_directories == "Y" or suppress_other_directories == "y":
                for i in kmer:
                    if i != kmer_choosen:  # will be transformed in not in (TODO)
                        os.rmdir("kmer_{}".format(i))
        else :
            print("No kmer had a mean of phase 0 superior to the threshold")

    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))


distribution_phase_plot(26, "R3", ".", 100)
distribution_phase_plot(27, "R3", ".", 100)
distribution_phase_plot(28, "R3", ".", 100)
# main()