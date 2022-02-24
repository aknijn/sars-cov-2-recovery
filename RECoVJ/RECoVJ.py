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
import json
import mysql.connector
from mysql.connector import errorcode

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def isNewLineage(inLineage):
    isNew = False
    if not '*' in inLineage:
        config = configparser.ConfigParser()
        config.read(TOOL_DIR + '/../recovery.conf')
        dbhost = config['db']['host']
        dbdatabase = config['db']['database']
        dbuser = config['db']['user']
        dbpassword = config['db']['password']
        config = {
            'user': dbuser,
            'password': dbpassword,
            'host': dbhost,
            'database': dbdatabase
        }
        sql = ("select * from v_sarscov2_lineages where Lineages = '" + inLineage + "'")
        try:
            cnx = mysql.connector.connect(**config)
            cursor = cnx.cursor(buffered=True)
            cursor.execute(sql)
            result = cursor.fetchone()
            if result == None:
                isNew = True
            cursor.close()
        except mysql.connector.Error as err:
            print(err)
        else:
            cnx.close()
    return isNew

def checkNewMutation(inSpike):
    isNew = False
    New = ""
    if not inSpike == '=':
        config = configparser.ConfigParser()
        config.read(TOOL_DIR + '/../recovery.conf')
        dbhost = config['db']['host']
        dbdatabase = config['db']['database']
        dbuser = config['db']['user']
        dbpassword = config['db']['password']
        config = {
            'user': dbuser,
            'password': dbpassword,
            'host': dbhost,
            'database': dbdatabase
        }
        try:
            cnx = mysql.connector.connect(**config)
            cursor = cnx.cursor(buffered=True)
            result = cursor.callproc('p_NewMutations', (inSpike, ''))
            New = result[1]
            if New is not None:
                if ";" in New:
                    isNew = True
            cursor.close()
        except mysql.connector.Error as err:
            print(err)
        else:
            cnx.close()
    return isNew, New

def isNotificaVariant(inLineage, inSpike):
    isNotifica = False
    with open(TOOL_DIR + '/Variants-EW', 'r') as f:
        variants = f.read().splitlines()
    for variant in variants:
        if '+' in variant:
            lineage_spike = variant.split('+')
            # AY.5.2 (== AY.5.2 or * or AY.*)
            if ((lineage_spike[0] == inLineage) or (lineage_spike[0] == '*') or ('*' in lineage_spike[0] and lineage_spike[0].replace('*','') in inLineage)) and (lineage_spike[1] in inSpike):
                isNotifica = True
                break
        else:
            if ((variant == inLineage) or ('*' in variant and variant.replace('*','') in inLineage)):
                isNotifica = True
                break
    return isNotifica

def getVariant(inLineage, inClade, inSpike, inLibrary):
    typeVariant = '-'
    if inLibrary != 'sang':
        typeVariant = getVariant_Lineage_Clade(inLineage, inSpike, 'Lineages')
    if typeVariant == '-':
        typeVariant = getVariant_Lineage_Clade(inClade, inSpike, 'Clades')
    return typeVariant

def getVariant_Lineage_Clade(inLineage_Clade, inSpike, inType):
    outVariant = '-'
    with open(TOOL_DIR + '/Variants-' + inType, 'r') as f:
        lines = f.read().splitlines()
    for line in lines:
        relations = line.split('\t')
        if relations[1] == '*':
            if ((relations[0] == inLineage_Clade) or ('*' in relations[0] and relations[0].replace('*','') in inLineage_Clade)):
                outVariant = relations[2]
                break
        else:
            if relations[0] == inLineage_Clade and (relations[1] in inSpike):
                outVariant = relations[2]
                break
    return outVariant

def colindex(gene):
    colindexes = {
        'orf1ab_nsp1' : 0,
        'orf1ab_nsp2' : 0,
        'orf1ab_nsp3' : 0,
        'orf1ab_nsp4' : 0,
        'orf1ab_nsp5' : 0,
        'orf1ab_nsp6' : 0,
        'orf1ab_nsp7' : 0,
        'orf1ab_nsp8' : 0,
        'orf1ab_nsp9' : 0,
        'orf1ab_nsp10' : 0,
        'orf1ab_nsp11' : 0,
        'orf1ab_nsp12' : 0,
        'orf1ab_nsp13' : 0,
        'orf1ab_nsp14' : 0,
        'orf1ab_nsp15' : 0,
        'orf1ab_nsp16' : 0,
        'S' : 1,
        'ORF3a' : 2, 
        'E' : 3,
        'M' : 4,
        'ORF6' : 5,
        'ORF7a' : 6,
        'ORF7b' : 7,
        'ORF8' : 8,
        'N' : 9,
        'ORF10' : 10
    }
    return colindexes.get(gene, 11)

def format_variants(str_variants):
    if str_variants == '=' or str_variants == 'ND':
        formatted_variants = str_variants
    else:
        formatted_variants = str_variants[1:-1]
    return formatted_variants
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strain', dest='strain', help='strain name')
    parser.add_argument('--librarytype', dest='librarytype', help='library type')
    parser.add_argument('--region', dest='region', help='region')
    parser.add_argument('--year', dest='year', help='year')
    parser.add_argument('--lineage', dest='lineage', help='pangolin')
    parser.add_argument('--clade', dest='clade', help='nextclade')
    parser.add_argument('--variants', dest='variants', help='Spike muations')
    parser.add_argument('--consensus', dest='consensus', help='consensus')
    parser.add_argument('--recovery_json', dest='recovery_json', help='output json')
    
    args = parser.parse_args()
    try:
        report_data = {}
        report_variants = []
        # prepare JSON output file
        report_data["information_name"] = args.strain
        report_data["region"] = args.region
        report_data["year"] = args.year
        # library type
        library = open(args.librarytype).readline().rstrip()
        if library == 'iont':
            report_data["sequence"] = "Ion Torrent"
        elif library == 'illu':
            report_data["sequence"] = "Illumina"
        elif library == 'nano':
            report_data["sequence"] = "Nanopore"
        elif library == 'sang':
            report_data["sequence"] = "Sanger"
        elif library == 'cons':
            report_data["sequence"] = "Consensus"
        # Ns in consensus
        with open(args.consensus, 'r') as cons_in:
            temp = cons_in.read().splitlines()
            consensus="".join(temp[1:])
        percN = (100.0 * consensus.count('N')) / (len(consensus))
        report_data["N_consensus"] = str(consensus.count('N')) + " (" + "{:.1f}".format(percN) + "%)"
        # obtain lineage and quality control from pangolin result and from Ns in consensus
        with open(args.lineage, 'r') as table_in:
            tab_lineage = [[str(col).rstrip() for col in row.split(',')] for row in table_in]
        report_data["lineage"] = tab_lineage[1][1]
        lineage = tab_lineage[1][1]
        if tab_lineage[1][len(tab_lineage[1])-2] != 'passed_qc':
            report_data["qc_status"] = 'Failed'
        else:
            if percN > 5.0:
                report_data["qc_status"] = 'Failed'
            else:
                report_data["qc_status"] = 'Passed'
        # obtain clade from nextclade result
        with open(args.clade, 'r') as table_in:
            tab_clade = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
        report_data["clade"] = tab_clade[1][1].strip('\"')
        # variants
        with open(args.variants, 'r') as table_in:
            tab_variants = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
        if library == 'sang':
            strDefault = "ND"
        else:
            strDefault = "="
        for i in range(0,12):
            report_variants.append(strDefault)
        report_variants[1] = "="
        for variant in tab_variants:
            if len(variant)>1:
                if variant[1] != 'Position' and colindex(variant[0]) != 11 and variant[6] != 'S':
                    if "DELETION" in variant[4]:
                        report_variants[colindex(variant[0])] = report_variants[colindex(variant[0])] + variant[6] + "_del; "
                    else:
                        report_variants[colindex(variant[0])] = report_variants[colindex(variant[0])] + variant[6] + "; "
        report_data["ORF1ab"] = format_variants(report_variants[0])
        report_data["S-protein"] = format_variants(report_variants[1])
        report_data["ORF3a"] = format_variants(report_variants[2])
        report_data["E-protein"] = format_variants(report_variants[3])
        report_data["M-protein"] = format_variants(report_variants[4])
        report_data["ORF6"] = format_variants(report_variants[5])
        report_data["ORF7a"] = format_variants(report_variants[6])
        report_data["ORF7b"] = format_variants(report_variants[7])
        report_data["ORF8"] = format_variants(report_variants[8])
        report_data["N-protein"] = format_variants(report_variants[9])
        report_data["ORF10"] = format_variants(report_variants[10])
        # Variante
        if isNewLineage(report_data["lineage"]):
            report_data["notifica"] = "nuovo lignaggio"
        else:
            isNewMutation, NewMutation = checkNewMutation(report_data["S-protein"])
            if isNewMutation:
                report_data["notifica"] = "nuova mutazione " + NewMutation
            else:
                report_data["notifica"] = "-"
        if isNotificaVariant(report_data["lineage"], report_data["S-protein"]):
            report_data["notifica"] = "Si"
        report_data["variante"] = getVariant(report_data["lineage"], report_data["clade"], report_data["S-protein"], library)
    finally:
        report = open(args.recovery_json, 'w')
        report.write("[" + json.dumps(report_data) + "]")
        report.close()

if __name__ == "__main__":
    main()
