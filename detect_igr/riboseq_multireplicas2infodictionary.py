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

def read_gff(file):
    dico = {}
    with open(file,'r') as f:
        for line in f:
            if line.startswith('chr') == False:
                continue
            
            ID        = line.split()[8].split(';')[0].split('=')[1]
            Parent    = line.split()[8].split(';')[1].split('=')[1]
            try:
                Overlap   = line.split()[8].split(';')[7].split('=')[1]
                Location  = line.split()[8].split(';')[5].split('=')[1]
            except:
                Overlap  = ''
                Location = ''
            dico[ID]  = {}
            dico[ID]['Parent']     = Parent
            dico[ID]['Overlap']    = Overlap
            dico[ID]['Location']   = Location
            dico[ID]['Intersect']  = ''
            dico[ID]['Replicas']   = {}
    return(dico)

def coverage_rate(coverage):
		""" Return the rate of position with at least a read in coverage (list).
		"""
		return round(sum((1 for x in coverage if x > 0)) / len(coverage), 5)
        
try:
    replicas_path  =    sys.argv[sys.argv.index("-replicas_path")+1]
    #replicas_path  = '/home/c.papadopoulos/Documents/data/riboseq/28mers/recalc/All_bams'
    gff_file       =    sys.argv[sys.argv.index("-gff")+1]
    #gff_file       = '/home/c.papadopoulos/Documents/data/gff_files/IGORF_noRNA_snoRNA_snRNA_transposable.gff'
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

print('Initialize the information dictionary')
info_dico   =   read_gff(gff_file)
print('Infromation dictionary created')
print('\n')
print('Read the GFF file')
GFF =  Gff(gff_file, all_as_high=True)
print('GFF file read \t DONE')


replicas = sorted(os.listdir(replicas_path))
r_count = 0
for replica in replicas:
    if replica.endswith('.bam'):
        r_count += 1
        replica_name = replica.split('.')[0]
        replica_code = 'R'+str(r_count)
        f_count = 0
        print('\n')
        #print('\t',replica_code)
        for feature in sorted(GFF.all_features(cast_into="Igorf")):
            f_count += 1
            print('\r\t' + str(replica_code) + '\t:\t' + str(f_count)+'\tsequences\' coverage', end = '')
            info_dico[feature.ID]['Replicas'][replica_code] = {}
            info_dico[feature.ID]['Replicas'][replica_code]['replica_name'] = replica_name
            
            info_dico[feature.ID]['Replicas'][replica_code]['All_frame'] = coverage_rate(feature.pos_cov(replicas_path+'/'+replica)['coverages'])
            # Number of total reads for this IGORF in rate_coverage
            info_dico[feature.ID]['Replicas'][replica_code]['Reads'] = sum((x for x in feature.pos_cov(replicas_path+'/'+replica)['coverages']))
            # Calcul of coverage by frame This value is not used. I use also count raw number
            coverage_by_frame = feature.frames_coverage(replicas_path+'/'+replica)
            info_dico[feature.ID]['Replicas'][replica_code]['0'] = coverage_rate(coverage_by_frame[0])
            info_dico[feature.ID]['Replicas'][replica_code]['1'] = coverage_rate(coverage_by_frame[1])
            info_dico[feature.ID]['Replicas'][replica_code]['2'] = coverage_rate(coverage_by_frame[2])
            
            info_dico[feature.ID]['Replicas'][replica_code]['positions_0'] = coverage_by_frame[0]
            info_dico[feature.ID]['Replicas'][replica_code]['positions_1'] = coverage_by_frame[1]
            info_dico[feature.ID]['Replicas'][replica_code]['positions_2'] = coverage_by_frame[2]

            # Number of total reads allowed value used
            # at frame 0
            info_dico[feature.ID]['Replicas'][replica_code]['covframe0'] = sum(coverage_by_frame[0])
            # at frame 1
            info_dico[feature.ID]['Replicas'][replica_code]['covframe1'] = sum(coverage_by_frame[1])
            # at frame 2
            info_dico[feature.ID]['Replicas'][replica_code]['covframe2'] = sum(coverage_by_frame[2])

    pickle.dump(info_dico, open(output_path + replica_code + '.fapick','wb'))

for ID in info_dico:
    reads_counter = 0
    for rep_code in info_dico[ID]['Replicas']:
        if info_dico[ID]['Replicas'][rep_code]['covframe0'] > 0:
            reads_counter += 1
    info_dico[ID]['Intersect'] = reads_counter
        



my_date = datetime.datetime.now()
out_name = output_name + '_' + str(datetime.datetime.now()).split()[0] + '.fapick'
pickle.dump(info_dico, open(output_path + out_name,'wb'))                                   








