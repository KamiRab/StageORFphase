#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:55:12 2019

@author: christospapadopoulos
"""

import sys,os
import pickle
import matplotlib
matplotlib.interactive(False)
import matplotlib.pyplot as plt

ribodico_file  = sys.argv[sys.argv.index("-ribodico")+1]
#ribodico_file = '/Volumes/LACIE SHARE/Phd/CDS_new_view/CDS_analysis_pipeline/Scer_transcriptome_all_replicas_2019-07-23.fapick'
ribo_dico = pickle.load(open( ribodico_file, "rb" ))

rep_counts = {}
for gene in ribo_dico:
    for replica in ribo_dico[gene]['Replicas']:
        
        if replica not in rep_counts.keys():
            rep_counts[replica] = {}
        
            rep_counts[replica]['reads_P0'] = 0
            rep_counts[replica]['reads_P1'] = 0
            rep_counts[replica]['reads_P2'] = 0
        
        rep_counts[replica]['reads_P0'] = rep_counts[replica]['reads_P0'] + int(ribo_dico[gene]['Replicas'][replica]['covframe0'])
        rep_counts[replica]['reads_P1'] = rep_counts[replica]['reads_P1'] + int(ribo_dico[gene]['Replicas'][replica]['covframe1'])
        rep_counts[replica]['reads_P2'] = rep_counts[replica]['reads_P2'] + int(ribo_dico[gene]['Replicas'][replica]['covframe2'])

perc_reads_in_phases = {}
for replica in rep_counts:
    perc_reads_in_phases[replica] = {}
    perc_reads_in_phases[replica]['P0'] = rep_counts[replica]['reads_P0'] / (rep_counts[replica]['reads_P0'] + rep_counts[replica]['reads_P1'] + rep_counts[replica]['reads_P2'])
    perc_reads_in_phases[replica]['P1'] = rep_counts[replica]['reads_P1'] / (rep_counts[replica]['reads_P0'] + rep_counts[replica]['reads_P1'] + rep_counts[replica]['reads_P2'])
    perc_reads_in_phases[replica]['P2'] = rep_counts[replica]['reads_P2'] / (rep_counts[replica]['reads_P0'] + rep_counts[replica]['reads_P1'] + rep_counts[replica]['reads_P2'])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(x=['P0','P1','P2'],height=[perc_reads_in_phases[replica]['P0'],perc_reads_in_phases[replica]['P1'],perc_reads_in_phases[replica]['P2']])
    ax.set_title(replica)
    ax.set_ylim(0,1)
    plt.show()
    fig.savefig(replica+'.png')
    plt.close("all")


