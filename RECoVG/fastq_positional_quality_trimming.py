# -*- coding: utf-8 -*-
"""
This tool trims paired-end FASTQ files on the basis of quality score or left/right position, retaining mate integrity.
Reads without mate after filtering are saved in a separate output file.
"""

import math
import optparse
import sys

def average(values):
    """ Arithmetic mean of a list of values """
    return math.fsum(values) / len(values) if len(values) else float('nan')


def phred2sanger(phred_scores):
    """ Convert an array of Phred quality scores (integers in [0, 93]) to a Sanger-encoded quality string"""
    return ''.join([chr(score + 33) for score in phred_scores])


def sanger2phred(sanger_string):
    """ Convert a Sanger-encoded quality string to an array of Phred quality scores"""
    return [ord(ch) - 33 for ch in sanger_string]


def trimming(sequ, qual, maxlengthtrim, lefttrim, righttrim, minqualtrim, avgqualtrim):
    """ Trimming of sequence and quality of a read"""
    # Maximum length trimming
    if maxlengthtrim != -1:
        if len(sequ) > maxlengthtrim:
            sequ = sequ[: maxlengthtrim]
            qual = qual[: maxlengthtrim]
    # Left- and right-side trimming
    if righttrim == 0:
        sequ = sequ[lefttrim :]
        qual = qual[lefttrim :]
    else:
        sequ = sequ[lefttrim : -righttrim]
        qual = qual[lefttrim : -righttrim]
    # Minimum quality right-side trimming
    while len(sequ) and qual[-1] < minqualtrim:
        qual = qual[:-1]
        sequ = sequ[:-1]
    # Average quality right-side trimming
    while len(sequ) and average(qual) < avgqualtrim:
        qual = qual[:-1]
        sequ = sequ[:-1]
    return sequ, qual


def __main__():
    parser = optparse.OptionParser()
    parser.add_option('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_option('-2', '--input2', dest='input2', help='reverse reads file in Sanger FASTQ format')
    parser.add_option('--maxlt', dest='maxlengthtrim', type='int', default=400, help='maximum length trimming')
    parser.add_option('--lt', dest='lefttrim', type='int', default=0, help='left-side trimming')
    parser.add_option('--rt', dest='righttrim', type='int', default=0, help='right-side trimming')
    parser.add_option('--minqt', dest='minqualtrim', type='int', default=15, help='minimum quality right-side trimming')
    parser.add_option('--avgqt', dest='avgqualtrim', type='float', default=20, help='average quality right-side trimming')
    parser.add_option('--minlf', dest='minlen', type='int', default=25, help='minimum length filtering')
    parser.add_option('--trimmed1', dest='trimmed1', help='trimmed forward FASTQ file')
    parser.add_option('--trimmed2', dest='trimmed2', help='trimmed reverse FASTQ file')
    parser.add_option('--trimmedunpaired', dest='trimmedunpaired', help='trimmed unpaired FASTQ file')
    parser.add_option('--log', dest='logfile', help='log file')
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')

    maxlengthtrim = options.maxlengthtrim
    lefttrim = options.lefttrim
    righttrim = options.righttrim
    minqualtrim = options.minqualtrim
    avgqualtrim = options.avgqualtrim
    minlen = options.minlen

    total_reads = 0
    discarded_reads = 0
    forward = open(options.input1)
    if options.input2:
        paired = True
        passing_paired_reads = 0
        passing_unpaired_reads = 0
        reverse = open(options.input2)
        trimmed_reverse = open(options.trimmed2, 'w')
        trimmed_unpaired = open(options.trimmedunpaired, 'w')
    else:
        paired = False
        passing_reads = 0
    trimmed_forward = open(options.trimmed1, 'w')
    log = open(options.logfile, 'w')
    try:
        while True:
            headL = next(forward).rstrip()
            sequL = next(forward).rstrip()
            commL = next(forward).rstrip()
            sangL = next(forward).rstrip()
            qualL = sanger2phred(sangL)
            trimmed_sequL, trimmed_qualL = trimming(sequL, qualL, maxlengthtrim, lefttrim, righttrim, minqualtrim, avgqualtrim)
            if paired:
                try:
                    headR = next(reverse).rstrip()
                    sequR = next(reverse).rstrip()
                    commR = next(reverse).rstrip()
                    sangR = next(reverse).rstrip()
                except StopIteration:
                    sys.exit('Reverse FASTQ file contain less reads than forward FASTQ file.')
                qualR = sanger2phred(sangR)
                trimmed_sequR, trimmed_qualR = trimming(sequR, qualR, maxlengthtrim, lefttrim, righttrim, minqualtrim, avgqualtrim)
            # Filter by residual length
            if paired:
                if len(trimmed_sequL) >= minlen and len(trimmed_sequR) >= minlen:
                    trimmed_forward.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + phred2sanger(trimmed_qualL) + '\n')
                    trimmed_reverse.write(headR + '\n' + trimmed_sequR + '\n' + commR + '\n' + phred2sanger(trimmed_qualR) + '\n')
                    passing_paired_reads += 1
                elif len(trimmed_sequL) >= minlen and len(trimmed_sequR) < minlen:
                    trimmed_unpaired.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + phred2sanger(trimmed_qualL) + '\n')
                    passing_unpaired_reads += 1
                elif len(trimmed_sequL) < minlen and len(trimmed_sequR) >= minlen:
                    trimmed_unpaired.write(headR + '\n' + trimmed_sequR + '\n' + commR + '\n' + phred2sanger(trimmed_qualR) + '\n')
                    passing_unpaired_reads += 1
                else:
                    discarded_reads += 1
            else:
                if len(trimmed_sequL) >= minlen:
                    trimmed_forward.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + phred2sanger(trimmed_qualL) + '\n')
                    passing_reads += 1
                else:
                    discarded_reads += 1
            total_reads += 1
    except StopIteration:
        if paired:
            try:
                next(reverse)
            except StopIteration:
                log.write("Total paired reads     : %d\n" % total_reads)
                log.write("Passing paired reads   : %d\n" % passing_paired_reads)
                log.write("Passing unpaired reads : %d\n" % passing_unpaired_reads)
                log.write("Discarded paired reads : %d\n" % discarded_reads)
            else:
                sys.exit('Forward FASTQ file contain less reads than reverse FASTQ file.')
        else:
            log.write("Total reads     : %d\n" % total_reads)
            log.write("Passing reads   : %d\n" % passing_reads)
            log.write("Discarded reads : %d\n" % discarded_reads)
    finally:
        forward.close()
        trimmed_forward.close()
        log.close()
        if paired:
            reverse.close()
            trimmed_reverse.close()
            trimmed_unpaired.close()


if __name__ == "__main__":
    __main__()
