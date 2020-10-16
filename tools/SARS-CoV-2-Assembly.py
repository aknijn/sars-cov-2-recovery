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
import sys
import os
import subprocess

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', dest='library', help='library type')
    parser.add_argument('-1', '--input1', dest='input1', help='trimmed forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('-2', '--input2', dest='input2', help='trimmed reverse reads file in Sanger FASTQ format')
    parser.add_argument('--spades_contigs', dest='spades_contigs', help='SPAdes contigs file')
    parser.add_argument('--spades_contig_stats', dest='spades_contig_stats', help='SPAdes contig stats file')
    parser.add_argument('--spades_scaffolds', dest='spades_scaffolds', help='SPAdes scaffolds file')
    parser.add_argument('--spades_scaffold_stats', dest='spades_scaffold_stats', help='SPAdes scaffold stats file')
    parser.add_argument('--spades_log', dest='spades_log', help='SPAdes log file')
    args = parser.parse_args()

    # Read library type from file and assemble input accordingly
    with open(args.library) as f:
        librarytype = f.readline().strip()
    # ION Torrent
    if librarytype=='iont':
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl " + args.spades_contigs + " " + args.spades_contig_stats + " " + args.spades_scaffolds + " " + args.spades_scaffold_stats + " " + args.spades_log + " NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} --iontorrent -s fastq:" + args.input1, shell=True)
    # Illumina
    elif librarytype=='illu':
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl " + args.spades_contigs + " " + args.spades_contig_stats + " " + args.spades_scaffolds + " " + args.spades_scaffold_stats + " " + args.spades_log + " NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} --pe1-ff --pe1-1 fastq:" + args.input1 + " --pe1-2 fastq:" + args.input2, shell=True)
    # Nanopore
    elif librarytype=='nano':
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl " + args.spades_contigs + " " + args.spades_contig_stats + " " + args.spades_scaffolds + " " + args.spades_scaffold_stats + " " + args.spades_log + " NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} -k '21,33,55,77,99,127' -s fastq:" + args.input1, shell=True)

if __name__ == "__main__":
    __main__()
