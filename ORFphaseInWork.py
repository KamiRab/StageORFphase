import argparse
import concurrent.futures.process
import errno
import os
import subprocess
import time
import pandas as pd
from datetime import datetime


def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORF phase')
    parser.add_argument("-gff",
                        type=str,
                        required=True,
                        # nargs="*",
                        help="GFF annotation file")
    parser.add_argument("-fasta",
                        type=str,
                        required=True,
                        # nargs="*",
                        help="FASTA file containing the genome sequence")
    parser.add_argument("-kmer",
                        type=list,
                        required=False,
                        nargs=2,
                        default=[26, 30])
    parser.add_argument("-riboseq",
                        type=str,
                        required=True,
                        # nargs="*",
                        help="FASTQ file containing the riboseq reads")
    parser.add_argument("-cutdir",
                        type=str,
                        required=False,  # if directory not given, launch cutadapt
                        # nargs="?",
                        help="Directory containing the directories of cutadapt files"
                        )
    parser.add_argument("-adapt",
                        type=str,
                        required=False,  # if not given, will look in list and if not found error message
                        nargs="*",
                        help="Nucleotidic sequence of the adaptor for riboseq")
    parser.add_argument("-thr",
                        type=int,
                        required=False,
                        nargs="?",
                        default=90,
                        help="Threshold for keeping the reads in phase 0")
    args = parser.parse_args()
    return args


def classify_reads(kmer, fastq_file, adaptor, rname,cutdir):
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
            -m {} -M {} -e 0.12 -o {}/kmer_{}/{}_kmer{}.fastq".format(fastq_file, adaptor, str(kmer), str(kmer),cutdir,
                                                                     str(kmer),
                                                                     rname, str(kmer))
    print("We are cutting the reads in {}-mers".format(kmer))
    # process =subprocess.run(cutadapt_cmd, shell=True)
    process = subprocess.run(cutadapt_cmd, shell=True, capture_output=True, universal_newlines=True)
    print(process.stdout, process.stderr)


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
    # a. Map the reads on the transcriptome and generate a sam file
    cmd_bowtie = "bowtie --wrapper basic-0 -v 2 -y -a -m 1 --best --strata -p 18 -x {}_CDS {}/kmer_{}/{}_kmer{}.fastq -S kmer_{}/{}_kmer{}.sam".format(
        gname, cutdir, str(kmer), rname, str(kmer), str(kmer), rname, str(kmer))
    # b. Transform sam into bam
    cmd_sam_bam = "samtools view -h -b -S {}/kmer_{}/{}_kmer{}.sam > {}/kmer_{}/{}_kmer{}.bam".format(cutdir, str(kmer), rname,
                                                                                                str(kmer),cutdir,
                                                                                                str(kmer), rname,
                                                                                                str(kmer))
    # c. Sort the bam file
    cmd_sam_sort = "samtools sort {}/kmer_{}/{}_kmer{}.bam -o {}/kmer_{}/{}_kmer{}_sorted.bam".format(cutdir, str(kmer), rname,
                                                                                                str(kmer),cutdir,
                                                                                                str(kmer), rname,
                                                                                                str(kmer))
    # d. Keep only the mapped reads
    cmd_sam_map = "samtools view -F 4 -f 0,16 -h -b {}/kmer_{}/{}_kmer{}_sorted.bam > {}/kmer_{}/{}_kmer{}_sorted_mapped.bam".format(
        cutdir, str(kmer), rname, str(kmer), cutdir, str(kmer), rname, str(kmer))
    # e. Index the bam file
    cmd_sam_index = "samtools index {}/kmer_{}/{}_kmer{}_sorted_mapped.bam".format(cutdir, str(kmer), rname, str(kmer))

    # generate count table
    cmd_count = "python ./detect_igr/BAM2Reads.py -bam kmer_{}/{}_kmer{}_sorted_mapped.bam -gff {}_transcriptome.gff " \
                "-output_path kmer_{} -output_name {}_kmer_{}".format(
        str(kmer), str(rname), str(kmer), str(gname), str(kmer), str(rname), str(kmer))

    process1 = subprocess.run(cmd_bowtie, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_bowtie, process1.stdout, process1.stderr)

    process2 = subprocess.run(cmd_sam_bam, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_sam_bam, process2.stdout, process2.stderr)

    process3 = subprocess.run(cmd_sam_sort, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_sam_sort, process3.stdout, process3.stderr)

    process4 = subprocess.run(cmd_sam_map, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_sam_map, process4.stdout, process4.stderr)

    process5 = subprocess.run(cmd_sam_index, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_sam_index, process5.stdout, process5.stderr)

    process6 = subprocess.run(cmd_count, shell=True, capture_output=True, universal_newlines=True)
    print("Command launched: ", cmd_count, process6.stdout, process6.stderr)

    print("The count table has been generated for kmer {}".format(kmer))


def perc_phase(kmer, rname, cutdir):
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
    phase0sum = tab["Number of p0"].sum()
    phase1sum = tab["Number of p1"].sum()
    phase2sum = tab["Number of p2"].sum()
    sumreads = phase0sum + phase1sum + phase2sum
    perc.append(phase0sum / sumreads)
    perc.append(phase1sum / sumreads)
    perc.append(phase2sum / sumreads)
    print(
        "For the {}_kmer_{}, there are : \n{:.3f}% reads in phase 0 \n{:.3f}% reads in phase 1 \n{:.3f}% reads in phase 2".format(
            rname, kmer, perc[0], perc[1], perc[2]))
    return perc


parameters = get_args()
kmer = range(parameters.kmer[0], parameters.kmer[1] + 1)
fastq_file = parameters.riboseq
genome_file = parameters.fasta
gff_file = parameters.gff

adaptor = parameters.adapt
rname = os.path.basename(fastq_file)

rname = os.path.splitext(rname)[0]
gname = os.path.basename(genome_file)
gname = os.path.splitext(gname)[0]
orfget_cmd = "orfget -fna {} -gff {} -features_include CDS -o {}_CDS -type nucl".format(
    genome_file, gff_file, gname)

print("Extract the Transcriptome:")
print("Launch :|t", orfget_cmd)
os.system(orfget_cmd)

transcriptome = read_multiFASTA(fasta_file=gname + "_CDS.nfasta")
with open(gname + "_transcriptome.gff", "w") as wgff:
    for i in transcriptome:
        wgff.write('{:20s}\t{}\t{:20s}\t{:d}\t{:d}\t{}\t{}\t{}\t{}\n'.format(
            i, 'SGD', 'gene', 1, len(transcriptome[i]), '.', '+', '0', str('ID=' + i)))

try:
    import cutadapt

    print("cutadadapt\t:\tCHECK")
except ImportError:
    print("cutadapt is not installed, please check")  # add link to installation
with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
    executor.map(classify_reads, kmer, [fastq_file] * len(kmer), [adaptor[0]] * len(kmer), [rname] * len(kmer))

# 3. We do the mapping of the reads (of 5 different sizes) on the CDS NT sequences
print("We start the mapping")
# a. Build the index:
cmd_A = 'bowtie-build {}_CDS.nfasta {}_CDS'.format(gname, gname)
process = subprocess.Popen(cmd_A, stderr=None, shell=True, stdout=subprocess.PIPE)
# os.system(cmd_A)
#
# while not os.path.exists(gname + "_CDS.1.ebwt"):
#     time.sleep(1)

####### check if bowtie and samtools are installed ?
# with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
#     executor.map(map_in_bam_and_count, kmer, [gname] * len(kmer), [rname] * len(kmer))

p0_by_kmer = {}
best_kmer = []
for i in kmer:
    p0_by_kmer[i] = perc_phase(i, rname, parameters.cutdir)
    if p0_by_kmer[i][0] > parameters.thr / 100:
        best_kmer.append(i)
print(best_kmer)


def main():
    start_time = datetime.now()
    global parameters
    parameters = get_args()
    kmer = range(parameters.kmer[0], parameters.kmer[1] + 1)
    fastq_file = parameters.riboseq
    genome_file = parameters.fasta
    gff_file = parameters.gff
    adaptor = parameters.adaptor

    rname = os.path.basename(fastq_file)
    rname = os.path.splitext(rname)[0]

    gname = os.path.basename(genome_file)
    gname = os.path.splitext(gname)[0]

    cutadapt_dir = "."
    # 1. ORFGET
    orfget_cmd = "orfget -fna {} -gff {} -features_include CDS -o {}_CDS -type nucl".format(
        genome_file, gff_file, gname)

    transcriptome = read_multiFASTA(fasta_file=gname + "_CDS.nfasta")
    with open(gname + "_transcriptome.gff", "w") as wgff:
        for i in transcriptome:
            wgff.write('{:20s}\t{}\t{:20s}\t{:d}\t{:d}\t{}\t{}\t{}\t{}\n'.format(
                i, 'SGD', 'gene', 1, len(transcriptome[i]), '.', '+', '0', str('ID=' + i)))

    print("Extract the Transcriptome:")
    print("Launch :|t", orfget_cmd)
    os.system(orfget_cmd)
    if not parameters.cutdir and not parameters.riboseq:
        print("You need to either give a fastq file or the directory containing the cutadapt files")
    elif not parameters.cutdir:
        try:
            import cutadapt
            print("The cutadapt process is launching. The reads are cut in {}-kmers".format(parameters.kmer))
        except ImportError:
            print('''cutadapt is not installed, please check ''')  # add link for installation
        with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
            executor.map(classify_reads, kmer, [fastq_file] * len(kmer), [adaptor[0]] * len(kmer),
                         [rname] * len(kmer))
    else:
        cutadapt_dir = parameters.cutdir
        for i in kmer:
            if not os.path.isdir("{}/kmer_{}".format(parameters.dir, i)):
                print("The directories need to be named as kmer_x with x the size of the kmer.")
                exit()

    print("We start the mapping")

    cmd_A = 'bowtie-build {}_CDS.nfasta {}_CDS'.format(gname, gname)
    os.system(cmd_A)
    while not os.path.exists(gname + "_CDS.1.ebwt"):
        time.sleep(1)

    with concurrent.futures.process.ProcessPoolExecutor(max_workers=None) as executor:
        executor.map(map_in_bam_and_count, kmer, [gname] * len(kmer), [rname] * len(kmer))
    end_time = datetime.now()
    print("\n\n")
    print('Duration: {}'.format(end_time - start_time))
