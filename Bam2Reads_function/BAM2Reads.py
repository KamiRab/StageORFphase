#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:04:13 2019

@author: c.papadopoulos
"""

# todo add headers
import argparse
import concurrent.futures.process
import os, sys
import re

from .loaders import Gff
import pickle
from datetime import datetime


def get_args():
    """

    :return: parameters
    """
    parser = argparse.ArgumentParser(description='Bam2Read')
    parser.add_argument("-gff",
                        type=str,
                        required=True,
                        nargs=1,
                        help="GFF annotation file of your regions of interest")
    parser.add_argument("-bam",
                        type=str,
                        required=True,
                        nargs=1,
                        help="FASTA file containing the genome sequence")
    parser.add_argument("-features_include",
                        type=str,
                        required=False,
                        nargs="*",
                        default=['all'],
                        help="Annotation features to be considered (By definition is all)")
    parser.add_argument("-features_exclude",
                        type=str,
                        required=False,
                        nargs="*",
                        default=["None"],
                        help="Annotation features not to be considered (By definition is None)")
    parser.add_argument("-outpath",
                        type=str,
                        # action='store',
                        required=True,
                        nargs=1,
                        help="The path to save the dictionary output")
    parser.add_argument("-outname",
                        type=str,
                        # action='store',
                        required=True,
                        nargs=1,
                        help="The name you want the output to have")
    args = parser.parse_args()
    return args


def coverage_rate(coverage):
    """ Return the rate of position with at least a read in coverage (list).
    """
    return round(sum((1 for x in coverage if x > 0)) / len(coverage), 5)


def count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, bam):
    with open(file_output + "_reads.tab", "w") as wtab, \
            open(file_output + "_periodicity_start.tab", "w") as wpstart, \
            open(file_output + "_periodicity_stop.tab", "w") as wpstop, \
            open(file_output+ "_periodicity_all.tab","w") as wper:
        start = datetime.now()
        features_found = 0
        wtab.write("ID\tNumber reads\tNumber p0\tNumber p1\tNumber p2\t"
                   "Perc. p0\tPerc. p1\tPerc. p2\n")
        wper.write("ID\tNumber p0\tNumber p1\tNumber p2\n")
        for x, feature in enumerate(gff_iterator):
            # Filtering of the features
            if re.search(elements_out, feature.ftype) and not re.search(elements_in, feature.ftype):
                continue

            if re.search(elements_out, feature.ftype) and re.search(elements_in, feature.ftype):
                print("{} include and exclude at the same time!".format(feature.ftype))
                exit()

            if not re.search(elements_in, feature.ftype) and elements_in != "(all)":
                continue
            features_found += 1
            # Coverage of the feature in the 3 frames
            coverage_by_frame = feature.frames_coverage(bam)
            reads_p0 = coverage_by_frame[0]
            reads_p1 = coverage_by_frame[1]
            reads_p2 = coverage_by_frame[2]

            # Number of reads by frame
            nb_reads_p0 = sum(coverage_by_frame[0])
            nb_reads_p1 = sum(coverage_by_frame[1])
            nb_reads_p2 = sum(coverage_by_frame[2])

            nb_reads_gene = nb_reads_p0 + nb_reads_p1 + nb_reads_p2

            # Percentage of reads by frame
            try:
                perc_reads_p0 = round(nb_reads_p0 / nb_reads_gene * 100, 2)
                perc_reads_p1 = round(nb_reads_p1 / nb_reads_gene * 100, 2)
                perc_reads_p2 = round(nb_reads_p2 / nb_reads_gene * 100, 2)
            except ZeroDivisionError:
                perc_reads_p0 = 0.0
                perc_reads_p1 = 0.0
                perc_reads_p2 = 0.0

            # Write on the reads count tab the total number of reads, the number and the percentage per phase for the feature
            wtab.write("{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{}\t{}\t{}\n".format(feature.ID, nb_reads_gene, nb_reads_p0,
                                                                           nb_reads_p1, nb_reads_p2, perc_reads_p0,
                                                                           perc_reads_p1, perc_reads_p2))

            #Periodicity of the feature of all length
            for aa in range(len(reads_p0)):
                wper.write("{}\t{}\t{}\t{}\t{}\n".format(feature.ID, aa+1, reads_p0[aa], reads_p1[aa], reads_p2[aa]))

            # Periodicity of the feature of length superior to 50
            if len(reads_p0) > 50:
                # We write the periodicity of the first 50 AA positions
                wpstart.write('{}\tp0\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p0[i]))
                wpstart.write('\n')
                wpstart.write('{}\tp1\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p1[i]))
                wpstart.write('\n')
                wpstart.write('{}\tp2\t'.format(feature.ID))
                for i in range(0, 51):
                    wpstart.write('{}\t'.format(reads_p2[i]))
                wpstart.write('\n')
                # We write the periodicity of the last 50 AA positions
                wpstop.write('{}\tp0\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p0[i]))
                wpstop.write('\n')
                wpstop.write('{}\tp1\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p1[i]))
                wpstop.write('\n')
                wpstop.write('{}\tp2\t'.format(feature.ID))
                for i in range(len(reads_p0) - 50, len(reads_p0)):
                    wpstop.write('{}\t'.format(reads_p2[i]))
                wpstop.write('\n')
        end = datetime.now()
    return "Duration : {}\n {} features corresponding the selection.\nThe mean time per selected feature is : {}\n".format(end-start, features_found,(end-start)/features_found)


def BAM2Reads(rname, gff_file, kmer, outname, elements_in=None, elements_out=None):
    if elements_in is None:
        elements_in = ["all"]
    if elements_out is None:
        elements_out = ["None"]

    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    print('Read the GFF file')
    start_time_gff = datetime.now()
    gff = Gff(gff_file, all_as_high=True)
    gff_iterator = sorted(gff.all_features(cast_into="Igorf"))
    end_time_gff = datetime.now()
    output_log = ""
    output_log = output_log+ 'Duration gff: {}'.format(end_time_gff - start_time_gff)
    print('GFF file read \t DONE')
    for size in kmer:
        if outname == "phasing":  # ORFphase
            file_input = "./kmer_{}/{}_kmer_{}_phasing".format(str(size), rname, str(size))
        else:  # ORFribomap
            file_input = "./kmer_{}/{}_kmer_{}_all".format(str(size), rname, str(size))
        bam_file = file_input + "_sorted_mapped.bam"
        file_output = "./kmer_{}/{}_kmer_{}_{}".format(str(size), rname, str(size), outname)
        output_log = output_log + "Kmer : {}\n".format(size)+count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, bam_file)
    end_time_all = datetime.now()
    return output_log+ '\n\nDuration b2r: {}'.format(end_time_all - start_time_gff)


def main():
    parameters = get_args()
    print('Read the GFF file')
    GFF = Gff(parameters.gff, all_as_high=True)
    gff_iterator = sorted(GFF.all_features(cast_into="Igorf"))
    print('GFF file read \t DONE')
    elements_in = "(" + ")|(".join(parameters.features_include) + ")"
    elements_out = "(" + ")|(".join(parameters.features_exclude) + ")"
    file_output = parameters.outpath + "/" + parameters.outname

    return count_percentage_reads_to_file(file_output, elements_in, elements_out, gff_iterator, parameters.bam)


if __name__ == "__main__":
    main()