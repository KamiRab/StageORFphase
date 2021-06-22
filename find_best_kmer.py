import argparse
import os

import pandas as pd
def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='ORFphase')
    parser.add_argument("-fastq",
                        type=str,
                        required=True,
                        nargs="?",
                        help="Riboseq file whose best phasing size is to find")
    parser.add_argument("-mode",
                        type=str,
                        choices=["best","thr"],
                        required=False,
                        default="best",
                        nargs="?",
                        help='Give kmer with best kmer ("best") or all the kmer over the threshold("thr")')
    parser.add_argument("-calc",
                        type=str,
                        choices=["median","mean"],
                        required=False,
                        default="mean",
                        nargs="?",
                        help='Find size by comparing mean or median')
    parser.add_argument("-thr",  # can only choose 1 threshold for all the files
                        type=int,
                        required=False,
                        nargs=1,
                        default=80,
                        help="Threshold for the percentage of good phased reads to have in percentage")
    args = parser.parse_args()
    return args
parameters = get_args()
riboseq_name = parameters.fastq
if riboseq_name.endswith(".fastq"):
    riboseq_name = os.path.basename(riboseq_name)
    riboseq_name = riboseq_name.split(".")[0]
mode = parameters.mode
kmer_tab = pd.read_table("{}_phase_median_mean.tab".format(riboseq_name))
if mode == "thr":
    thr = parameters.thr/100
    if parameters.calc == "mean":
        size = kmer_tab[kmer_tab["Mean_p0"]>thr]
        print(kmer_tab[kmer_tab["Mean_p0"]>thr].Size.values)
    else:
        print(kmer_tab[kmer_tab["Median_p0"]>thr].Size.values)
else:
    if parameters.calc == "mean":
        print(kmer_tab.iloc[kmer_tab.Mean_p0.argmax(),1])
    else:
        print(kmer_tab.iloc[kmer_tab.Median_p0.argmax(),1])
