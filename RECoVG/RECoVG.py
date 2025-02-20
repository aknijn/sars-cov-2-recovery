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
import shutil
import subprocess

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--library', dest='library', help='library type')
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('-2', '--input2', dest='input2', help='reverse reads file in Sanger FASTQ format')
    parser.add_argument('--covidref_aligned', dest='covidref_aligned', help='COVIDRef aligned sorted bam file')
    parser.add_argument('--reference_genbank', dest='reference_genbank', help='reference genbank file')
    parser.add_argument('--reference_fasta', dest='reference_fasta', help='reference fasta file')
    parser.add_argument('--proteinentcovid19', dest='proteinentcovid19', help='ProteineNt_Covid19 fasta file')
    args = parser.parse_args()

    os.system("ln -s " + os.popen("which trimmomatic.jar").read().strip() + " trimmomatic.jar")
    # Elaborate each library type appropriately
    # ION Torrent
    if args.library=='iont':
        # TRIMMING
        subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 '" + args.input1 + "' trimmed1.fq LEADING:15 TRAILING:15 MINLEN:30 SLIDINGWINDOW:4:15", shell=True)
        # FILTER HUMAN GENOME
        subprocess.call("bowtie2 -p ${GALAXY_SLOTS:-4} -x '" + TOOL_DIR + "/data/Humangenome' -U trimmed1.fq --un filtered.fq --very-fast | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o human_genome_aligned", shell=True)
        # ALIGN SARS-COV-2 GENOME
        subprocess.call("bowtie2 -p ${GALAXY_SLOTS:-4} -x '" + TOOL_DIR + "/data/genome' -U filtered.fq --very-sensitive-local | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o " + args.covidref_aligned, shell=True)
    # Illumina
    elif args.library=='illu':
        # TRIMMING
        # subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar PE -threads ${GALAXY_SLOTS:-6} -phred33 '" + args.input1 + "' '" + args.input2 + "' trimmed1.fq fastq_out_r1_unpaired trimmed2.fq fastq_out_r2_unpaired LEADING:15 TRAILING:15 MINLEN:30 SLIDINGWINDOW:4:15", shell=True)
        subprocess.call("python " + TOOL_DIR + "/fastq_positional_quality_trimming.py --maxlt -1 --lt 0 --rt 5 --minqt 15 --avgqt 15.0 --minlf 30 -1 '" + args.input1 + "' --trimmed1 trimmed1.fq --log log.txt -2 '" + args.input2 + "' --trimmed2 trimmed2.fq --trimmedunpaired trimmedunpaired.fq", shell=True)
        # FILTER HUMAN GENOME
        subprocess.call("bowtie2 -p ${GALAXY_SLOTS:-4} -x '" + TOOL_DIR + "/data/Humangenome' -1 trimmed1.fq -2 trimmed2.fq --un-conc filtered.fq --very-fast | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o human_genome_aligned", shell=True)
        # ALIGN SARS-COV-2 GENOME
        subprocess.call("bowtie2 -p ${GALAXY_SLOTS:-4} -x '" + TOOL_DIR + "/data/genome' -1 filtered.1.fq -2 filtered.2.fq --very-sensitive-local | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o aligned.bam", shell=True)
        # TRIM PRIMERS
        subprocess.call("ivar trim -e -m 30 -q 20 -s 4 -i aligned.bam -b " + TOOL_DIR + "/data/QIAseq_artic.bed -p primer_trimmed", shell=True)
        subprocess.call("samtools sort -@ \${GALAXY_SLOTS:-1} -o primer_trimmed.sorted.bam primer_trimmed.bam", shell=True)
        shutil.copy("primer_trimmed.sorted.bam", args.covidref_aligned)
    # Nanopore
    elif args.library=='nano':
        if args.input1.lower().endswith(".fast5"):
            os.system("ln -s " + args.input1 + " chopped.fq")
        else:
            # PORECHOP
            subprocess.call("porechop -i '" + args.input1 + "' --format 'fastq' -o 'chopped.fq'", shell=True)
        # TRIMMING
        subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 '" + args.input1 + "' trimmed1.fq LEADING:9 TRAILING:9 MINLEN:30", shell=True)
        # FILTER HUMAN GENOME
        subprocess.call("bowtie2 -p ${GALAXY_SLOTS:-4} -x '" + TOOL_DIR + "/data/Humangenome' -U trimmed1.fq --un filtered.fq --very-fast | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o human_genome_aligned", shell=True)
        # ALIGN SARS-COV-2 GENOME
        subprocess.call("minimap2 -t ${GALAXY_SLOTS:-4} " + TOOL_DIR + "/data/genome.fa filtered.fq -a | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o " + args.covidref_aligned, shell=True)
    # Sanger
    elif args.library=='sang':
        # ALIGN SARS-COV-2 GENOME
        subprocess.call("minimap2 -t ${GALAXY_SLOTS:-4} " + TOOL_DIR + "/data/genome_Spike.fa '" + args.input1 + "' -a | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o " + args.covidref_aligned, shell=True)
    # Consensus
    elif args.library=='cons':
        # ALIGN SARS-COV-2 GENOME
        subprocess.call("sed '/^[^>]/s/-//g' '" + args.input1 + "' > no-fasta", shell=True)
        subprocess.call("sed '/^[^>]/s/[^ATGCatgc]/N/g' no-fasta > clean_fasta", shell=True)
        subprocess.call("minimap2 -t ${GALAXY_SLOTS:-4} " + TOOL_DIR + "/data/genome.fa clean_fasta -a | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o " + args.covidref_aligned, shell=True)

    # COPY CORRESPONDING REFERENCE (Deprecated)
    if args.library=='sang':
        shutil.copy(TOOL_DIR + "/data/CovidREF_Spike.gbk", args.reference_genbank)
        shutil.copy(TOOL_DIR + "/data/genome_Spike.fa", args.reference_fasta)
    else:
        shutil.copy(TOOL_DIR + "/data/CovidREF.gbk", args.reference_genbank)
        shutil.copy(TOOL_DIR + "/data/genome.fa", args.reference_fasta)
    shutil.copy(TOOL_DIR + "/data/ProteineNt_Covid19.fasta", args.proteinentcovid19)

if __name__ == "__main__":
    __main__()
