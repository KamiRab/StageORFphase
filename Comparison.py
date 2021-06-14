# gene en phase !=
# par rapport à pourcentage normale ?
# somme réplicat ?
import os

import pandas as pd
def phase_decision_median(cutdir, kmer, riboseq_name, thr):
    p0_by_kmer = {}
    kmer_over_thr = {}
    best_size = {0: [0, 0], 1: [0, 0], 2: [0, 0]}
    median={}
    for size in kmer:
        median[size]= []
        tab = pd.read_table("{}/kmer_{}/{}_kmer_{}_phasing_reads.tab".format(cutdir, size, riboseq_name, size),
                            sep='\t',
                            names=["ID", "Number of reads", "Number of p0", "Number of p1", "Number of p2",
                                   "Percentage of p0", "Percentage of p1", "Percentage of p2"])
        # Sum of phase columns in tab
        phase0sum = tab["Number of p0"].sum()
        phase1sum = tab["Number of p1"].sum()
        phase2sum = tab["Number of p2"].sum()
        sumreads = tab["Number of reads"].sum()

        # List with percentage of each phase
        median[size].append(phase0sum / sumreads)
        median[size].append(phase1sum / sumreads)
        median[size].append(phase2sum / sumreads)
        # print(
        #     "For the {}_kmer_{}, there are : \n{:.3f}% reads in phase 0 \n{:.3f}% reads in phase 1 \n"
        #     "{:.3f}% reads in phase 2".format(riboseq_name, size, p0_by_kmer[size][0],
        #                                       p0_by_kmer[size][1], p0_by_kmer[size][2]))
        if median[size][0] > thr / 100:
            kmer_over_thr[size] = median[size][0]
        for phase in range(3):
            if median[size][phase] > best_size[phase][1]:
                best_size[phase] = [size, median[size][phase]]
    print("For {} the kmers with a mediaan superior to {}% are: ".format(riboseq_name, thr), kmer_over_thr)
    print("The best size for {} for each phase is:".format(riboseq_name), best_size)
    if not kmer_over_thr:
        print("No kmer had a mediaan of phase 0 superior to the threshold for {}".format(riboseq_name))
    return [kmer_over_thr, best_size]
def main():
    kmer = range(26,31)
    riboseq_tab = pd.read_csv("PolyAribo.csv",sep=",", names = ["GSE","SRR", "adapt"])
    for SRR in  riboseq_tab["SRR"]:
        GSE = riboseq_tab[riboseq_tab["SRR"]==SRR].GSE.values[0]
        for size in kmer:
            best_median_over_thr = phase_decision_median(".", kmer, SRR, 80)[0]
            best_median_by_phase = phase_decision_median(".", kmer, SRR, 80)[1]
            if not os.path.isfile("best_median_size_polyA.tab".format(GSE)):
                with open("best_median_size_polyA_{}.tab".format(GSE), "w") as median_tab:
                    median_tab.write("GSE\tRiboseq\tPhase\tSize\tMedian\n")
            else:
                with open("best_median_size_polyA_{}.tab".format(GSE), "a") as median_tab:
                    for phase in range(3):
                        median_tab.write(
                            "{}\t{}\t{}\t{}\t{}\n".format(GSE, SRR, phase,
                                                          best_median_by_phase[phase][0], best_median_by_phase[phase][1]))

main()