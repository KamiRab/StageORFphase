import errno
import os
import subprocess
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
    if len(gff) == len(fasta) and len(fasta) == len(fastq):
        if parameters.adapt and not parameters.cutdir:
            adapt = parameters.adapt
            if not len(fastq) == len(adapt):
                if len(adapt) == 1:
                    print("You have given only one adaptor. All the fastq are considered as having the same adaptor")
                else:
                    print("You need to put the same number of adapters and files")
                    exit()
            print("The associations of the files are:")
            for i in range(len(gff)):
                print("GFF: {} \tFASTA: {}\tFASTAQ: {}\t adaptor: {}"
                      .format(gff[i], fasta[i], fastq[i], adapt[i]))
        elif not parameters.adapt and parameters.cutdir:
            cutdir = parameters.cutdir
            if not len(fastq) == len(cutdir):
                if len(cutdir) == 1:
                    print("You have given only one cutadapt directory. All the cutadapt results have to be inside")
                else:
                    print("You need to put the same number of cutadapt directories as files")
                    exit()
            print("The associations of the files are:")
            for i in range(len(gff)):
                print("GFF: {} \tFASTA: {}\tFASTAQ: {}\t cutadapt directory: {}"
                      .format(gff[i], fasta[i], fastq[i], cutdir[i]))
        elif not parameters.adapt and not parameters.cutdir:
            print("You need to enter either the cutadapt directory or the adaptor sequence(s)")
            exit()
    else:
        print("There is something wrong with the number of files input. There must be the same number of gff, "
              "fasta files and fastq files")
        exit()

def cut_reads(kmer, fastq_file, adaptor, rname, cutdir):
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
        os.mkdir("{}/kmer_{}".format(cutdir, str(kmer)))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    cutadapt_cmd = "cutadapt {} -j 1 --quality-base=33 -a {} \
-m {} -M {} -e 0.12 -o {}/kmer_{}/{}_kmer{}.fastq".format(fastq_file, adaptor, str(kmer), str(kmer), cutdir,
                                                          str(kmer),
                                                          rname, str(kmer))
    print("We are cutting the reads in {}-mers".format(kmer))
    process = subprocess.run(cutadapt_cmd, shell=True)


def map_in_bam_and_count(kmer, gname, rname, cutdir, output_name, gff_to_compare):
    '''
    Map the differents sizes of reads on the genome and
    :param kmer:
    :param gname:
    :param rname:
    :param cutdir:
    :return:
    '''

    file_output = "{}/kmer_{}/{}_kmer{}".format(cutdir, str(kmer), rname, str(kmer))
    ebwt_basename = gname + "_" + output_name
    # Commands
    # b. Map the reads on the Bowtie index and generate a sam file of the found alignments
    cmd_bowtie = "bowtie --wrapper basic-0 -v 2 -y -a -m 1 --best --strata -p 18 -x {} {}.fastq " \
                 "-S {}.sam".format(ebwt_basename, file_output, file_output)
    # c. Transform sam into bam file
    cmd_sam_bam = "samtools view -h -b -S {}.sam > {}.bam".format(file_output, file_output)
    # d. Sort the bam file
    cmd_sam_sort = "samtools sort {}.bam -o {}_sorted.bam".format(file_output, file_output)
    # e. Keep only the mapped reads
    cmd_sam_map = "samtools view -F 4 -f 0,16 -h -b {}_sorted.bam " \
                  "> {}_sorted_mapped.bam".format(file_output, file_output)
    # f. Index the bam file
    cmd_sam_index = "samtools index {}_sorted_mapped.bam".format(file_output)

    # g. generate count table
    cmd_count = "python ./detect_igr/BAM2Reads.py -bam {}_sorted_mapped.bam -gff " \
                "{} -output_path {}/kmer_{} -output_name {}_kmer_{}".format(file_output, gff_to_compare,
                                                                            cutdir, str(kmer), rname, str(kmer))

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

def reads_phase_plot(table, kmer, rname, reads_thr):
    '''

    :param table:
    :param kmer:
    :param rname:
    :param reads_thr:
    :return:
    '''
    tab_select = table[table['Number of reads'] > reads_thr]
    plt.figure()
    tab_select.boxplot(column=["Percentage of p0", "Percentage of p1", "Percentage of p2"])
    # plt.show()
    plt.savefig('./Boxplot_phases_{}_CDS_kmer{}.png'.format(rname, kmer))


def reads_periodicity(kmer, rname, cutdir, type):
    '''

    :param kmer:
    :param rname:
    :param cutdir:
    :return:
    '''
    tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_periodicity_{}.tab".format(cutdir, kmer, rname, kmer, type), sep='\t',
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
    sns_plot.savefig("{}/kmer_{}/{}_kmer_{}_periodicity_{}.png".format(cutdir, kmer, rname, kmer, type))
