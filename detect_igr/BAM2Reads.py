#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:04:13 2019

@author: c.papadopoulos
"""


# todo add parser orfget
import argparse
import concurrent.futures.process
import os,sys
import re

from loaders import Gff
import pickle
import datetime
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
                        #action='store',
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
        
# try:
#     BAM            =    sys.argv[sys.argv.index("-bam")+1]
#     gff_file       =    sys.argv[sys.argv.index("-gff")+1]
#     output_path    =    sys.argv[sys.argv.index("-output_path")+1]
#     output_name    =    sys.argv[sys.argv.index("-output_name")+1]
# except:
#     print('''
#     This is my usefull help!!!
#
#     You need to give some inputs first:
#         -replicas_path : The path where ALL the .bam files of riboseq are stored
#         -gff           : The GFF annotation file of your regions of interest
#         -output_path   : The path to save the dictionary output
#         -output_name   : The name you want your output to have!!!
#                          ATTENTION!!!
#                          The name will not be Exectly as you give it:
#                              it will have also the date of generation & will
#                              look like:
#                                  NAME_yyyy_mm_dd.fapick
#     ''')
#     exit()





def count_percentage_reads_to_file(file_output, elements_in, elements_out, gff, bam):
    with open(file_output + "_reads.tab", "w") as wtab, open(
            file_output + "_periodicity_start.tab", "w") as wpstart, open(
            file_output + "_periodicity_stop.tab", "w") as wpstop:
        for x, feature in enumerate(sorted(gff.all_features(cast_into="Igorf"))):
            # f_count += 1
            # print('\r\t' + str(replica_code) + '\t:\t' + str(f_count)+'\tsequences\' coverage', end = '')
            if re.match(elements_out, feature.ftype) and not re.match(elements_in, feature.ftype):
                continue

            if re.match(elements_out, feature.ftype) and re.match(elements_in, feature.ftype):
                print("{} include and exclude at the same time!".format(feature.ftype))
                exit()

            if not re.match(elements_in, feature.ftype) and elements_in != "(all)":
                continue
            coverage_by_frame = feature.frames_coverage(bam)
            reads_p0 = coverage_by_frame[0]
            reads_p1 = coverage_by_frame[1]
            reads_p2 = coverage_by_frame[2]

            nb_reads_p0 = sum(coverage_by_frame[0])
            nb_reads_p1 = sum(coverage_by_frame[1])
            nb_reads_p2 = sum(coverage_by_frame[2])

            nb_reads_gene = nb_reads_p0 + nb_reads_p1 + nb_reads_p2

            try:
                perc_reads_p0 = round(nb_reads_p0 / nb_reads_gene * 100, 2)
                perc_reads_p1 = round(nb_reads_p1 / nb_reads_gene * 100, 2)
                perc_reads_p2 = round(nb_reads_p2 / nb_reads_gene * 100, 2)
            except ZeroDivisionError:
                perc_reads_p0 = 0.0
                perc_reads_p1 = 0.0
                perc_reads_p2 = 0.0

            wtab.write("{:20s}\t{:d}\t{:d}\t{:d}\t{:d}\t{}\t{}\t{}\n".format(feature.ID, nb_reads_gene, nb_reads_p0,
                                                                             nb_reads_p1, nb_reads_p2, perc_reads_p0,
                                                                             perc_reads_p1, perc_reads_p2))

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

def BAM2Reads(cutdir, rname, gff_file, kmer, elements_in=None, elements_out=None):
    if elements_in is None:
        elements_in = ["all"]
    if elements_out is None:
        elements_out = ["None"]

    elements_in = "(" + ")|(".join(elements_in) + ")"
    elements_out = "(" + ")|(".join(elements_out) + ")"
    print('Read the GFF file')
    gff = Gff(gff_file, all_as_high=True)
    print('GFF file read \t DONE')
    with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
        for size in kmer:
            file_output = "{}/kmer_{}/{}_kmer{}".format(cutdir, str(size), rname, str(size))
            bam_file = file_output + "_sorted_mapped.bam"
            executor.submit(count_percentage_reads_to_file, file_output, elements_in, elements_out, gff, bam_file)

# BAM2Reads("..", "R3", "../mapping_orf_Saccer3.gff", [27,28])

parameters = get_args()
print('Read the GFF file')
GFF = Gff(parameters.gff, all_as_high=True)
print('GFF file read \t DONE')
elements_in = "(" + ")|(".join(parameters.features_include) + ")"
elements_out = "(" + ")|(".join(parameters.features_exclude) + ")"
file_output = parameters.outpath + "/" + parameters.outname

count_percentage_reads_to_file(file_output,elements_in,elements_out,GFF,parameters.bam)




