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
import configparser
import sys
import os
import shutil
import subprocess
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../RECoVLibs/")
from recovdb import IridaDb

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def getMetadata(inputfiles, inuser, inspecies):
    iridaDb = IridaDb(inspecies)
    with open(inputfiles, 'r') as f:
        content = f.readlines()
    idfiles = [getIdFile(x.rstrip('\n')) for x in content]
    files_id = ",".join(idfiles)
    records = iridaDb.metadata_for_consensus(inuser, files_id)
    iridaDb.close()
    return records

# Obtain idFile from file path
def getIdFile(filename):
    splitFilename = filename.split("/")
    if (splitFilename[5][0]=='A'):
         inIdFile = splitFilename[6]
    else:
         inIdFile = splitFilename[5]
    return inIdFile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', dest='input_files', help='input files')
    parser.add_argument('--user', dest='user', help='user')
    parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--multifasta', dest='multifasta', help='multifasta')
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(TOOL_DIR + '/../recovery.conf')
    iridadir = config['fs']['output_path']
    metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), args.species)
    if metadata:
        # create MultiFasta
        subprocess.run("touch multifasta", shell=True)
        for meta_row in metadata:
            with open(iridadir + "/" + meta_row[0], 'r') as cons_in:
                temp = cons_in.read().splitlines()
                consensus="".join(temp[1:])
            with open("multifasta", 'a') as cons_out:
                cons_out.write(">" + meta_row[1] +"\n")
                cons_out.write(consensus +"\n")
        shutil.copyfile("multifasta", args.multifasta)
    else:
        with open(args.multifasta, 'w') as multi_fasta:
            multi_fasta.write("Errore\n")

if __name__ == "__main__":
    main()


