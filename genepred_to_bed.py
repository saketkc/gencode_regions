#!/usr/bin/env python
"""
Extract first/last exon and CDS from GenePred formatted file
"""

from __future__ import print_function
import pandas
import argparse
import sys

def extract_first_coding_exon(row):
    cdsStart = int(row['cdsStart'])
    cdsEnd = int(row['cdsEnd'])
    strand = row['strand']
    exonStarts = row['exonStarts']
    exonEnds = row['exonEnds']

    ## Noncoding
    if cdsStart == cdsEnd:
        return

    exonStarts_all = exonStarts.split(',')
    exonEnds_all = exonEnds.split(',')
    exonEnds_all = [int(x)  for x in exonEnds_all if x]
    exonStarts_all = [int(x) for x in exonStarts_all if x]
    name = '{}'.format(row['name'])
    if strand =='+':
        for i in range(0, len(exonEnds_all)):
            ## If the end of exon is greater than the start of CDS
            if exonEnds_all[i]>=cdsStart:
                return '{}\t{}\t{}\t{}\t.\t{}'.format(row['chrom'], cdsStart, exonEnds_all[i], name, strand)
    if strand == '-':
        for i in range(len(exonEnds_all)-1, 0, -1):
            if (exonStarts_all[i] <= cdsEnd):
                return '{}\t{}\t{}\t{}\t.\t{}'.format(row['chrom'], exonStarts_all[i], cdsEnd, name, strand)

def extract_last_coding_exon(row):
    cdsStart = int(row['cdsStart'])
    cdsEnd = int(row['cdsEnd'])
    strand = row['strand']
    exonStarts = row['exonStarts']
    exonEnds = row['exonEnds']

    ## Noncoding
    if cdsStart == cdsEnd:
        return

    exonStarts_all = exonStarts.split(',')
    exonEnds_all = exonEnds.split(',')
    exonEnds_all = [int(x)  for x in exonEnds_all if x]
    exonStarts_all = [int(x) for x in exonStarts_all if x]
    name = '{}'.format(row['name'])
    if strand =='+':
        for i in range(len(exonEnds_all)-1, 0, -1):
            ## If the end of exon is greater than the start of CDS
            if exonEnds_all[i]<=cdsEnd and exonStarts_all[i]>=cdsStart:
                return '{}\t{}\t{}\t{}\t.\t{}'.format(row['chrom'], exonStarts_all[i], exonEnds_all[i], name, strand)
    if strand == '-':
        for i in range(0, len(exonStarts_all)):
            if exonStarts_all[i]>=cdsStart and exonEnds_all[i]<=cdsEnd:
                return '{}\t{}\t{}\t{}\t.\t{}'.format(row['chrom'], exonStarts_all[i], exonEnds_all[i], name, strand)


def get_CDS(row):
    cdsStart = int(row['cdsStart'])
    cdsEnd = int(row['cdsEnd'])
    strand = row['strand']
    name = '{}'.format(row['name'])
    if (-cdsStart+cdsEnd)>0:
        return '{}\t{}\t{}\t{}\t.\t{}'.format(row['chrom'], cdsStart, cdsEnd, name, strand)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_exon', help='Output first exon BED', action='store_true')
    parser.add_argument('--last_exon', help='Output last exon BED', action='store_true')
    parser.add_argument('--cds', help='Output CDS BED', action='store_true')
    parser.add_argument('genepred', help='Path to GTF')

    args = parser.parse_args()
    if not (args.first_exon or args.last_exon or args.cds):
        sys.stderr.write('Should select one of --first_exon, --last_exon, --cds')
        sys.exit(1)

    if args.first_exon:
        fetch_func = extract_first_coding_exon
    elif args.last_exon:
        fetch_func = extract_last_coding_exon
    elif args.cds:
        fetch_func = get_CDS

    df = pandas.read_table(args.genepred, header=None)
    print(len(df.columns))
    if len(df.columns)==10:
        df.columns=['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
    elif len(df.columns)==15:
        df.columns=['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    else:
        raise RuntimeError('Input not in GenePred format')

    for index, row in df.iterrows():
        record = fetch_func(row)
        if record:
            print(record)

