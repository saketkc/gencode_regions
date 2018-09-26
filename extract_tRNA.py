#!/usr/bin/env
'''
Extract tRNA coordinates from GTF
'''
import sys
import GTF
import numpy as np
import pandas as pd

def main(GENCODE):
    gc = GTF.dataframe(GENCODE)
    gc.gene_id = gc.gene_id.replace(to_replace=r'\.[0-9]+', value='', regex=True)

    idx = (gc.feature == 'transcript') & gc.transcript_type.str.contains('tRNA')
    tRNA = gc.ix[idx, ['seqname','start','end','transcript_id','gene_name', 'strand']]
    tRNA.start = tRNA.start.astype(int)
    tRNA.end = tRNA.end.astype(int)
    tRNA.sort_values(by=['seqname','start','end'], inplace=True)
    tRNA.to_csv('tRNA_transcripts.bed', sep='\t', header=False, index=False)

    idx = (gc.feature == 'gene') & gc.gene_type.str.contains('tRNA')
    tRNA = gc.ix[idx, ['seqname','start','end','gene_id','gene_name', 'strand']]
    tRNA.start = tRNA.start.astype(int)
    tRNA.end = tRNA.end.astype(int)
    tRNA.sort_values(by=['seqname','start','end'], inplace=True)
    tRNA.to_csv('tRNA_genes.bed', sep='\t', header=False, index=False)

if __name__ == '__main__':
    main(sys.argv[1])

