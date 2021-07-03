#remove deletion/insertion caused by homopolymers and NGS errors
from Bio import SeqIO
import subprocess
import csv
import sys
import argparse

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_consensus', dest='first_consensus', help='first_consensus file')
    parser.add_argument('--reference_fasta', dest='reference_fasta', help='reference_fasta file')
    parser.add_argument('--minority_variants', dest='minority_variants', help='minority_variants file')
    parser.add_argument('--majority_variants', dest='majority_variants', help='majority_variants file')
    args = parser.parse_args()

    subprocess.call("cat " + args.reference_fasta + " " + args.first_consensus + " > sequences.fasta", shell=True)
    subprocess.call("mafft --quiet --auto sequences.fasta > all.fasta", shell=True)

    records=list(SeqIO.parse("all.fasta", "fasta"))
    reference=records[0].seq
    sequence=records[1].seq
    name_sequence=records[1].id

    csv_max_file = open(args.majority_variants)
    read_csv_max = list(csv.reader(csv_max_file, delimiter="\t"))
    csv_min_file= open(args.minority_variants)
    read_csv_min = list(csv.reader(csv_min_file, delimiter="\t"))
    
    read_csv_minmax=[]
    for line in read_csv_max[1:]:
        if line[5].find("FRAME_SHIFT")!=-1:
            read_csv_minmax.append(line)
    for line in read_csv_min[1:]:
        if line[5].find("FRAME_SHIFT")!=-1:
            read_csv_minmax.append(line)

    new_sequence=''
    i=0
    lunghezza=len(sequence)-1
    while i<lunghezza:
        if sequence[i-1]!='-' and sequence[i+1]!='-' and sequence[i]=='-':
            position_tab=i
            for line in read_csv_minmax:
                if position_tab==int(line[1]) and read_csv_minmax[read_csv_minmax.index(line)][1]!=read_csv_minmax[read_csv_minmax.index(line)-1][1]:
                    nucleotide=line[2][1].lower()
                    new_sequence+=nucleotide
            i+=1
        else:
            new_sequence+=sequence[i]
            i+=1
    
    i=0
    lunghezza=len(reference)
    to_remove=[]
    while i<lunghezza-1:
        if reference[i-1]!='-' and reference[i+1]!='-' and reference[i]=='-':
            to_remove.append(i)
            i+=1
        else:
            i+=1
    
    if len(to_remove)>=1:
        to_remove.sort(reverse=True)
        for i in to_remove:
            for line in read_csv_minmax:
                if i == int(line[1]):
                    new_sequence=new_sequence[:i]+new_sequence[i+1:]

    new_sequence=new_sequence.replace("-","")
    fasta=open("consensus.fasta", "w")
    fasta.write(">"+name_sequence+"\n")
    fasta.write(new_sequence.upper())
    fasta.close

if __name__ == "__main__":
    __main__()