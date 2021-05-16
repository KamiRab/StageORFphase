#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:04:13 2019

@author: c.papadopoulos
"""



import os,sys
from loaders import Gff
import pickle
import datetime

def coverage_rate(coverage):
		""" Return the rate of position with at least a read in coverage (list).
		"""
		return round(sum((1 for x in coverage if x > 0)) / len(coverage), 5)
        
try:
    BAM            =    sys.argv[sys.argv.index("-bam")+1]
    gff_file       =    sys.argv[sys.argv.index("-gff")+1]
    output_path    =    sys.argv[sys.argv.index("-output_path")+1]
    output_name    =    sys.argv[sys.argv.index("-output_name")+1]
except:
    print('''
    This is my usefull help!!!
    
    You need to give some inputs first:
        -replicas_path : The path where ALL the .bam files of riboseq are stored
        -gff           : The GFF annotation file of your regions of interest
        -output_path   : The path to save the dictionary output
        -output_name   : The name you want your output to have!!!
                         ATTENTION!!! 
                         The name will not be Exectly as you give it:
                             it will have also the date of generation & will
                             look like:
                                 NAME_yyyy_mm_dd.fapick
    ''')
    exit()


print('Read the GFF file')
GFF =  Gff(gff_file, all_as_high=True)
print('GFF file read \t DONE')



with open(output_path+'/'+output_name+"_reads.tab","w") as wtab, open(output_path+'/'+output_name+"_periodicity_start.tab","w") as wpstart,open(output_path+'/'+output_name+"_periodicity_stop.tab","w") as wpstop:
    for x,feature in enumerate(sorted(GFF.all_features(cast_into="Igorf"))):
        #f_count += 1
        #print('\r\t' + str(replica_code) + '\t:\t' + str(f_count)+'\tsequences\' coverage', end = '')
       
        coverage_by_frame = feature.frames_coverage(BAM)
        reads_p0 = coverage_by_frame[0]
        reads_p1 = coverage_by_frame[1]
        reads_p2 = coverage_by_frame[2]
        
        nb_reads_p0 = sum(coverage_by_frame[0])
        nb_reads_p1 = sum(coverage_by_frame[1])
        nb_reads_p2 = sum(coverage_by_frame[2])
        
        nb_reads_gene = nb_reads_p0 + nb_reads_p1 + nb_reads_p2
        
        try:
            perc_reads_p0 = round(nb_reads_p0/nb_reads_gene*100,2)
            perc_reads_p1 = round(nb_reads_p1/nb_reads_gene*100,2)
            perc_reads_p2 = round(nb_reads_p2/nb_reads_gene*100,2)
        except:
            perc_reads_p0 = 0.0
            perc_reads_p1 = 0.0
            perc_reads_p2 = 0.0
        
        wtab.write("{:20s}\t{:d}\t{:d}\t{:d}\t{:d}\t{}\t{}\t{}\n".format(feature.ID,nb_reads_gene,nb_reads_p0,nb_reads_p1,nb_reads_p2,perc_reads_p0,perc_reads_p1,perc_reads_p2))
        
        if len(reads_p0) > 50:
            # We write the periodicity of the first 50 AA positions
            wpstart.write('{}\tp0\t'.format(feature.ID))
            for i in range(0,51):
                wpstart.write('{}\t'.format(reads_p0[i]))
            wpstart.write('\n')
            wpstart.write('{}\tp1\t'.format(feature.ID))
            for i in range(0,51):
                wpstart.write('{}\t'.format(reads_p1[i]))
            wpstart.write('\n')
            wpstart.write('{}\tp2\t'.format(feature.ID))
            for i in range(0,51):
                wpstart.write('{}\t'.format(reads_p2[i]))
            wpstart.write('\n')
            # We write the periodicity of the last 50 AA positions
            wpstop.write('{}\tp0\t'.format(feature.ID))
            for i in range(len(reads_p0) - 50,len(reads_p0)):
                wpstop.write('{}\t'.format(reads_p0[i]))
            wpstop.write('\n')
            wpstop.write('{}\tp1\t'.format(feature.ID))
            for i in range(len(reads_p0) - 50,len(reads_p0)):
                wpstop.write('{}\t'.format(reads_p1[i]))
            wpstop.write('\n')
            wpstop.write('{}\tp2\t'.format(feature.ID))
            for i in range(len(reads_p0) - 50,len(reads_p0)):
                wpstop.write('{}\t'.format(reads_p2[i]))
            wpstop.write('\n')
            








