#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 09:54:00 2019

@author: christospapadopoulos
"""

import pickle,os,sys


def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return(unique_list)

def get_file_title(title):
    unique_elements = unique(title)
    final_title = ''
    for elem in unique_elements:
        final_title = final_title + elem + '\t'
    final_title = final_title + '\n'
    return(final_title)


infodico    =   sys.argv[sys.argv.index("-infodico")+1]
infotable   =   sys.argv[sys.argv.index("-infotable")+1]

try:
    info_dico = pickle.load(open( infodico, "rb" ))
    print(infodico,'\t read with success!')
except:
    print('Problem with reading infodico')
    exit()
    
title = ['ID']
with open( infotable , 'w' ) as fw:
    for i in info_dico : 
        fw.write('{:20s}\t'.format(i))
        point_one = int(i.split('_')[2].split('-')[0])
        point_two = int(i.split('_')[2].split('-')[1])
        start = min(point_one,point_two)
        stop  = max(point_one,point_two)
        width = abs(point_one - point_two)
        fw.write('{:7d}\t'.format(width))
        title.append('Size')
        fw.write('{:10d}\t{:10d}\t'.format(start,stop))
        title.append('Start')
        title.append('Stop')
        strand = i.split('_')[1]
        fw.write('{:3s}\t'.format(strand))
        title.append('Strand')
        fw.write('{:20s}\t'.format(info_dico[i]['Parent']))
        title.append('Parent')
        fw.write('{:10s}\t'.format(info_dico[i]['Location']))
        title.append('Location')
        fw.write('{:8s}\t'.format(info_dico[i]['Overlap']))
        title.append('Overlap')
        if info_dico[i]['Intersect'] == '':
            fw.write('{:5s}\t'.format('-'))
        else:
            fw.write('{:5d}\t'.format(info_dico[i]['Intersect']))
        title.append('Intersect')
        for rep in info_dico[i]['Replicas']:
            for opt in info_dico[i]['Replicas'][rep]:
                if opt == 'replica_name':
                    continue
                fw.write('{:10.5f}\t'.format(info_dico[i]['Replicas'][rep][opt]))
                title.append(str(opt)+'_'+str(rep))
        
        fw.write('\n')
            
        
        
final_title = get_file_title(title)         
with open( infotable , "r+") as f: s = f.read(); f.seek(0); f.write(final_title + s)
      
        
        