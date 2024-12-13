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
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../RECoVLibs/")
from recovdb import IridaDb

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def isNewLineage(inLineage):
    isNew = False
    if not '*' in inLineage:
        iridaDb = IridaDb('Coronavirus')
        isNew = iridaDb.isNewLineage(inLineage)
        iridaDb.close()
    return isNew

def checkNewMutation(inSpike):
    isNew = False
    New = ""
    if not inSpike == '=':
        iridaDb = IridaDb('Coronavirus')
        isNew, New = iridaDb.checkNewMutation(inSpike)
        iridaDb.close()
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

def isNotificaVariant2(inVariante, inSpike, inOrf1ab, inEprot):
    isNotifica = "-"
    if ((inVariante == 'Omicron') and ('I484V' in inOrf1ab) and ('A488S' in inOrf1ab)):
        isNotifica = "NSP2:I484V & NSP3:A488S"
    if ((inVariante == 'Omicron') and ('R346T' in inSpike) and ('K444' in inSpike) and ('N460K' in inSpike)):
         isNotifica = "S:R346T & S:K444 & S:N460K"
    if ((inVariante == 'Omicron') and ('E484K' in inSpike) and ('V445H' in inSpike)):
         isNotifica = "BA.2.86"
    if (('P191S' in inOrf1ab) and ('F225L' in inOrf1ab)):
         isNotifica = "ORF1a:P191S & ORF1a:F225L"
    if (('V642G' in inSpike) and ('S691P' in inSpike) and ('T791I' in inSpike)):
         isNotifica = "S:V642G & S:S691P & S:T791I"
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
        report_data["information_name"] = open(args.strain).readline().rstrip()
        report_data["region"] = open(args.region).readline().rstrip()
        report_data["year"] = open(args.year).readline().rstrip()
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
        consensusN = consensus.replace("n", "N")
        consensusN = consensusN.replace("?", "N")
        consensusN = consensusN.replace("!", "N")
        percN = (100.0 * consensusN.count('N')) / (len(consensusN))
        report_data["N_consensus"] = str(consensusN.count('N')) + " (" + "{:.1f}".format(percN) + "%)"
        # obtain lineage and quality control from pangolin result and from Ns in consensus
        with open(args.lineage, 'r') as table_in:
            tab_lineage = [[str(col).rstrip() for col in row.split(',')] for row in table_in]
        report_data["lineage"] = tab_lineage[1][1]
        # lineage = tab_lineage[1][1]
        if tab_lineage[1][len(tab_lineage[1])-3] != 'pass':
            report_data["qc_status"] = 'Failed'
        else:
            if percN > 5.0:
                report_data["qc_status"] = 'Failed'
            else:
                report_data["qc_status"] = 'Passed'
        # obtain clade from nextclade result
        with open(args.clade, 'r') as table_in:
            tab_clade = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
        report_data["clade"] = tab_clade[1][2].strip('\"')
        clade_lineage = tab_clade[1][3].strip('\"')
        if report_data["clade"] == "recombinant":
            report_data["clade"] = "recombinant " + clade_lineage
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
                if variant[1] != 'Position' and colindex(variant[0]) != 11 and len(variant) > 6:
                    if variant[len(variant)-1] != 'S':
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
        report_data["variante"] = getVariant(report_data["lineage"], report_data["clade"], report_data["S-protein"], library)
        # Notifica
        report_data["notifica"] = "-"
        if isNotificaVariant(report_data["lineage"], report_data["S-protein"]):
            report_data["notifica"] = "Si"
        #if report_data["variante"] != 'Omicron' and report_data["S-protein"].count(';') > 15:
        #    report_data["notifica"] = "Si"
        if report_data["lineage"][0:1] == "X":
            report_data["variante"] = "Ricombinante"
        if report_data["clade"][0:11] == "recombinant" and report_data["sequence"] != "Sanger" and report_data["ORF1ab"] != "=" and report_data["S-protein"] != "=" and report_data["qc_status"] != 'Failed':
            report_data["variante"] = "Ricombinante"
        elif report_data["clade"][0:11] == "recombinant" and report_data["sequence"] != "Sanger":
            report_data["notifica"] = "analisi incerta"
        elif report_data["clade"][0:11] == "recombinant":
            report_data["clade"] = "ND"
            report_data["notifica"] = "analisi incerta"
        notificaVariantSpot = isNotificaVariant2(report_data["variante"], report_data["S-protein"], report_data["ORF1ab"], report_data["E-protein"])
        if notificaVariantSpot != "-":
            report_data["notifica"] = notificaVariantSpot
        if isNewLineage(report_data["lineage"]):
            report_data["notifica"] = "nuovo lignaggio " + report_data["lineage"]
        else:
            isNewMutation, NewMutation = checkNewMutation(report_data["S-protein"])
            if isNewMutation:
                report_data["notifica"] = "nuova mutazione " + NewMutation
    finally:
        report = open(args.recovery_json, 'w')
        report.write("[" + json.dumps(report_data) + "]")
        report.close()

if __name__ == "__main__":
    main()
