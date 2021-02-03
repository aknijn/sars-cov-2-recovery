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
import json

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
    if str_variants == '=':
        formatted_variants = str_variants
    else:
        formatted_variants = str_variants[1:-1]
    return formatted_variants
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strain', dest='strain', help='strain')
    parser.add_argument('--region', dest='region', help='region')
    parser.add_argument('--year', dest='year', help='year')
    parser.add_argument('--lineage', dest='lineage', help='lineage')
    parser.add_argument('--variants', dest='variants', help='variants')
    parser.add_argument('--recovery_json', dest='recovery_json', help='recovery_json')
    
    args = parser.parse_args()
    try:
        report_data = {}
        report_variants = []
        # prepare JSON output file
        report_data["information_name"] = args.strain
        report_data["region"] = args.region
        report_data["year"] = args.year
        with open(args.lineage) as table_in:
            tab_lineage = [[str(col).rstrip() for col in row.split(',')] for row in table_in]
        report_data["lineage"] = tab_lineage[1][1] + " (" + tab_lineage[1][2] + ")"
        if tab_lineage[1][4] != 'passed_qc':
            report_data["qc_status"] = 'Failed'
        else:
            report_data["qc_status"] = 'Passed'
        with open(args.variants) as table_in:
            tab_variants = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]        
        for i in range(0,12):
            report_variants.append("=")
        for variant in tab_variants:
            if variant[1] != 'POS' and colindex(variant[0]) != 11:
                report_variants[colindex(variant[0])] = report_variants[colindex(variant[0])] + variant[7] + ";"
                # report_variants[colindex(variant[0])] = report_variants[colindex(variant[0])] + variant[1] + ":" + variant[2] + "/" + variant[3] + ";"
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
        # report_data["Intergenic"] = format_variants(report_variants[11])
    finally:
        report = open(args.recovery_json, 'w')
        report.write("[" + json.dumps(report_data) + "]")
        report.close()

if __name__ == "__main__":
    main()
