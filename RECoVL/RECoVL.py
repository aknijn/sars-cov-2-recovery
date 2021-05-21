#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import os
import sys
import subprocess
import csv

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def getAABase(AA):
    AABase=''
    for x in AA:
        if x.isdigit():
            AABase+=x
    return AABase

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants', dest='variants', help='spike variants')
    parser.add_argument('--strain', dest='strain', help='strain name')
    parser.add_argument('--lineage', dest='lineage', help='lineage')
    
    args = parser.parse_args()
    with open(TOOL_DIR + '/watchlist_unique.txt', 'r') as wl:
        mutations={}
        for l in wl:
            line=l.split(':')
            mutations[line[0]]=line[1].replace('\n', '')
    
    with open(args.variants, 'r') as spike:
        read_csv = list(csv.reader(spike, delimiter="\t"))
    spike_mut=[]
    for l in read_csv[1:]:
        if l[6].find(getAABase(l[6]))==1:
            spike_mut.append(l[6][1:])
        else:
            spike_mut.append(l[6])
    
    spike_mut_set=set(spike_mut)
    intersection={}
    for key, value in mutations.items():
        list2=set([s.replace(' ', '') for s in value.split(',')])
        int=list(spike_mut_set.intersection(list2))
        if len(int)>=1:
            intersection[key]=int

    intersection_sorted=sorted(intersection.items(), key=lambda x: len(x[1]), reverse=True)

    len_best_match = 0
    best_matches = []
    for pos in intersection_sorted:
        if len_best_match == 0:
            len_best_match = len(pos[1])
        if len(pos[1]) == len_best_match:
            best_matches.append(pos[0])
    if len(best_matches) == 0:
        best_matches.append("ND")
    elif len(best_matches) > 1 and "B.1" in best_matches:
        best_matches.remove("B.1")
    report = open(args.lineage, 'w')
    report.write("taxon,lineage,conflict,-,-,-,status,note\n")
    report.write(args.strain + ",(*)" + ";".join(best_matches) + ",ND (0),-,-,-,seq_len:,(*) lineage non può essere con certezza ottenuto da sequenza Sanger")
    report.close()

if __name__ == "__main__":
    main()
