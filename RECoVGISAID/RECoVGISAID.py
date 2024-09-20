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
import csv
from datetime import datetime

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.dirname(os.path.abspath(__file__)) + '/../RECoVLibs/'
CONFIG_FILE = TOOL_DIR + '/../recovery.conf'
if os.path.isdir('/mnt/pulsar/files')
    TOOL_DIR = '/mnt/pulsar/files'
    LIB_DIR = TOOL_DIR + '/RECoVLibs/'
    CONFIG_FILE = TOOL_DIR + '/recovery.conf'

sys.path.append(LIB_DIR)
from recovdb import IridaDb


def getMetadata(inputfiles, inuser, inspecies):
    iridaDb = IridaDb(inspecies)
    with open(inputfiles, 'r') as f:
        content = f.readlines()
    idfiles = [getIdFile(x.rstrip('\n')) for x in content]
    files_id = ",".join(idfiles)
    records = iridaDb.metadata_for_gisaid(inuser, files_id)
    header_file = iridaDb.gisaid_header_file
    iridaDb.close()
    with open(header_file) as csv_in:
        header = [str(col).rstrip() for col in csv_in.split(',')]
    return header, records

def updateGISAID(gisaid_file, sample_dict):
    iridaDb = IridaDb(inspecies)
    with open(gisaid_file, 'r') as f:
        uploaded_seqs = f.readlines()
    for uploaded_seq in uploaded_seqs:
        col_seq = uploaded_seq.split(': ')
        if len(col_seq) == 2:
            if col_seq[0][:8] == "SUCCESS;":
                iridaDb.update_externalId(str(col_seq[1]).rstrip(), sample_dict[col_seq[0][9:-22]])
    iridaDb.close()

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
    # parser.add_argument('--species', dest='species', help='species')
    parser.add_argument('--metadati', dest='metadati', help='metadati')
    parser.add_argument('--multifasta', dest='multifasta', help='multifasta')
    parser.add_argument('--response', dest='response', help='response')
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    iridadir = config['fs']['output_path']
    gisaid_username = config['gisaid']['username']
    gisaid_password = config['gisaid']['password']
    gisaid_clientid = config['gisaid']['clientid']

    species = "Coronavirus" #args.species
    if species == "Coronavirus" or species == "SARS-CoV-2":
        multifasta_prefix = 'icogen_'
        gisaid_cli = TOOL_DIR + '/covCLI'
    elif species == "Influenza A/Influenza B":
        multifasta_prefix = 'iflugen_'
        gisaid_cli = TOOL_DIR + '/fluCLI'
    else:
        multifasta_prefix = 'error_'
        gisaid_cli = 'error'

    csv_header, metadata = getMetadata(args.input_files, args.user.replace("__at__", "@"), species)
    sample_dict = {}
    if metadata:
        # Create MultiFasta
        multifasta = multifasta_prefix + datetime.today().strftime('%Y%m%d') + '.fa'
        subprocess.run("echo '' > " + multifasta, shell=True)
        for meta_row in metadata:
            sample_dict[meta_row[2]] = str(meta_row[-1])
            subprocess.run("echo '>" + meta_row[2] + "' >> " + multifasta, shell=True)
            subprocess.run("tail -n +2 " + iridadir + "/" + meta_row[0] + " >> " + multifasta, shell=True)
            subprocess.run("echo '' >> " + multifasta, shell=True)
        shutil.copyfile(multifasta, args.multifasta)
        # Create csv
        with open(args.metadati, 'w') as meta_csv:
            meta_csv_writer = csv.writer(meta_csv)
            meta_csv_writer.writerow(csv_header)
            for meta_row in metadata:
                meta_csv_writer.writerow(meta_row[1:-3])
        # Perform upload and create response file
        subprocess.run(gisaid_cli + " upload --username " + gisaid_username + " --password " + gisaid_password + " --clientid " + gisaid_clientid + " --fasta " + multifasta + " --metadata " + args.metadati + " --log " + args.response, shell=True)
        # Write GISAID codes to sample.externalId
        updateGISAID(args.response, sample_dict)
    else:
        with open(args.multifasta, 'w') as multi_fasta:
            multi_fasta.write("Errore\n")
        with open(args.metadati, 'w') as meta_csv:
            meta_csv.write("Errore\n")
        with open(args.response, 'w') as response_file:
            response_file.write("Errore\n")

if __name__ == "__main__":
    main()
