import errno
import os
import subprocess
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parameters_verification(parameters):
    '''
    Function to verify the number on input for fastq, fasta and gff are correct. And there is at least an adapter
    or cutadapt directory
    :param parameters:
    :return:
    '''
    gff = parameters.gff
    fasta = parameters.fasta
    fastq = parameters.fastq
    if len(gff) == len(fasta) and len(fasta) == len(fastq) or len(gff) == 1 or len(fasta) == 1:
        if len(gff) == 1 and len(fastq) > 1:
            print("Warning : you entered only one gff file, it will be used for all fastq file ")
            gff = gff * len(fastq)
        if len(fasta) == 1 and len(fastq) > 1:
            print("Warning : you entered only one fasta file, it will be used for all fastq file ")
            fasta = fasta * len(fastq)
        if parameters.adapt and not parameters.cutdir:
            adapt = parameters.adapt
            if not len(fastq) == len(adapt):
                if len(adapt) == 1:
                    print("You have given only one adaptor. All the fastq are considered as having the same adaptor")
                    adapt = adapt*len(fastq)
                else:
                    print("You need to put the same number of adapters and files")
                    exit()
            print("The associations of the files are:")
            for i in range(len(fastq)):
                print("GFF: {} \tFASTA: {}\tFASTAQ: {}\t adaptor: {}\n"
                      .format(gff[i], fasta[i], fastq[i], adapt[i]))
        elif not parameters.adapt and parameters.cutdir:
            cutdir = parameters.cutdir
            if not len(fastq) == len(cutdir):
                if len(cutdir) == 1:
                    print("You have given only one cutadapt directory. All the cutadapt results have to be inside")
                    cutdir = cutdir * len(fastq)
                else:
                    print("You need to put the same number of cutadapt directories as files")
                    exit()
            print("The associations of the files are:")
            for i in range(len(fastq)):
                print("GFF: {} \tFASTA: {}\tFASTAQ: {}\t cutadapt directory: {}\n"
                      .format(gff[i], fasta[i], fastq[i], cutdir[i]))
        elif not parameters.adapt and not parameters.cutdir:
            print("You need to enter either the cutadapt directory or the adaptor sequence(s)")
            exit()
    else:
        print("There is something wrong with the number of files input. There must be the same number of gff, "
              "fasta files and fastq files")
        exit()


def cut_reads(kmer, fastq_file, adaptor, rname):
    """
    Use cutadapt to trim the 3' adapter and filtering the reads to have only of one size
    :param kmer: filtering sizes
    :param fastq_file: file containings the riboseq reads
    :param adaptor: adaptor to trim
    :param rname: fastq_name without extension
    :param cutdir: directory that will contain the cutadapt output
    :return:
    """
    try:
        os.mkdir("./kmer_{}".format(str(kmer)))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    cutadapt_cmd = "cutadapt {} -j 1 --quality-base=33 -a {} \
-m {} -M {} -e 0.12 -o ./kmer_{}/{}_kmer_{}.fastq".format(fastq_file, adaptor, str(kmer), str(kmer),
                                                           str(kmer),
                                                           rname, str(kmer))
    # print("We are cutting the reads in {}-mers".format(kmer))
    process = subprocess.run(cutadapt_cmd, shell=True, universal_newlines=True, stdout=subprocess.PIPE)
    return process.stdout



def map2bam(kmer, gname, rname, output_name):
    '''
    Map the differents sizes of reads on the genome and creates bam file
    :param output_name:
    :param kmer:
    :param gname:
    :param rname:
    :param cutdir:
    :return:
    '''
    fastq_input = "./kmer_{}/{}_kmer_{}".format(str(kmer), rname, str(kmer))
    file_output = "./kmer_{}/{}_kmer_{}_{}".format(str(kmer), rname, str(kmer), output_name)
    ebwt_basename = gname + "_" + output_name
    # Commands
    # b. Map the reads on the Bowtie index and generate a sam file of the found alignments
    cmd_bowtie = "bowtie --wrapper basic-0 -v 2 -y -a -m 1 --best --strata -p 18 -x {} {}.fastq " \
                 "-S {}.sam".format(ebwt_basename, fastq_input, file_output)
    # c. Transform sam into bam file
    cmd_sam_bam = "samtools view -h -b -S {}.sam > {}.bam".format(file_output, file_output)
    # d. Sort the bam file
    cmd_sam_sort = "samtools sort {}.bam -o {}_sorted.bam".format(file_output, file_output)
    # e. Keep only the mapped reads
    cmd_sam_map = "samtools view -F 4 -f 0,16 -h -b {}_sorted.bam " \
                  "> {}_sorted_mapped.bam".format(file_output, file_output)
    # f. Index the bam file
    cmd_sam_index = "samtools index {}_sorted_mapped.bam".format(file_output)

    # Processes
    # print("Command launched: ", cmd_bowtie)
    process1 = subprocess.run(cmd_bowtie, shell=True, universal_newlines=True, stdout=subprocess.PIPE)

    print("Command launched: ", cmd_sam_bam)
    process2 = subprocess.run(cmd_sam_bam, shell=True, universal_newlines=True, stdout=subprocess.PIPE)

    print("Command launched: ", cmd_sam_sort)
    process3 = subprocess.run(cmd_sam_sort, shell=True, universal_newlines=True, stdout=subprocess.PIPE)

    print("Command launched: ", cmd_sam_map)
    process4 = subprocess.run(cmd_sam_map, shell=True, universal_newlines=True, stdout=subprocess.PIPE)

    print("Command launched: ", cmd_sam_index)
    process5 = subprocess.run(cmd_sam_index, shell=True, universal_newlines=True, stdout=subprocess.PIPE)
    return process1.stdout+process2.stdout+process3.stdout+process4.stdout+process5.stdout


def reads_phase_plot(table, kmer, rname, reads_thr, outname):
    '''

    :param outname:
    :param cutdir:
    :param table:
    :param kmer:
    :param rname:
    :param reads_thr:
    :return:
    '''
    tab_select = table[table['Number reads'] > reads_thr]
    plt.figure()
    tab_select.boxplot(column=["Perc. p0", "Perc. p1", "Perc. p2"])
    # plt.show()
    plt.savefig('./kmer_{}/{}_kmer_{}_{}_Boxplot_phases.png'.format(kmer, rname, kmer, outname))
    plt.figure()
    plt.title("Repartition of the phases for the reads of size {} ".format(kmer))
    sns.set_style('whitegrid')
    try:
        tab_select[["Perc. p0", "Perc. p1", "Perc. p2"]].plot.kde(bw_method=0.5)
    except ValueError as ve:
        print("The density plot can't be generated for kmer_{} :".format(kmer), ve)
    plt.savefig('./kmer_{}/{}_kmer_{}_{}_Density_phases.png'.format(kmer, rname, kmer, outname))
    plt.close('all')


def reads_periodicity(kmer, rname, outname, type):
    '''

    :param type:
    :param outname:
    :param kmer:
    :param rname:
    :param cutdir:
    :return:
    '''
    tab = pd.read_table("./kmer_{}/{}_kmer_{}_{}_periodicity_{}.tab".format(kmer, rname, kmer, outname, type), sep='\t',
                        header=None)
    tab.columns = [str(x) for x in range(len(tab.columns))]
    tab = tab.dropna(how="all", axis=1).drop(columns=["0"]).rename(columns={'1': "Phase"})
    tab_agg = pd.melt(tab, id_vars=["Phase"], var_name="Position", value_name="Number of reads")
    tab_agg["Position"] = pd.to_numeric(tab_agg["Position"]) - 2
    tab_agg = tab_agg.groupby(["Phase", "Position"], as_index=False).sum().sort_values(["Phase", "Position"])
    sns_plot = sns.catplot(x="Position", y="Number of reads", hue="Phase", data=tab_agg, kind="bar", height=8.27,
                           aspect=11.7 / 8.27)
    if type == "stop":
        sns_plot.set_xticklabels([x for x in range(-50, 0)])
    sns_plot.savefig("./kmer_{}/{}_kmer_{}_{}_periodicity_{}.png".format(kmer, rname, kmer, outname, type))
