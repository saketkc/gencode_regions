#!/usr/bin/env python
import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('seaborn-muted')
import pandas
import seaborn as sns
sns.set_style(style='white')

def load_df(filepath, name, cutoff):
    """Load bed file and return sites less
    than 'cutoff' length"""
    df = pandas.read_table(filepath, header=None, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    df['length'] = df.stop-df.start
    df['domain'] =  name
    df = df[df['length']<cutoff]
    return df

def plot_size_distribution(input_dir, title, outprefix):

    first_coding_exon = os.path.join(input_dir, 'first_exons.bed')
    last_coding_exon = os.path.join(input_dir, 'last_exons.bed')
    utr3 = os.path.join(input_dir, '3UTRs.bed')
    utr5 = os.path.join(input_dir, '5UTRs.bed')
    intron = os.path.join(input_dir, 'introns.bed')
    cds = os.path.join(input_dir, 'cds.bed')

    first_coding_exon_df = load_df(first_coding_exon, 'first_coding_exon', 1500)
    last_coding_exon_df = load_df(last_coding_exon, 'last_coding_exon', 1500)
    utr3_df = load_df(utr3, 'utr3', 1500)
    utr5_df = load_df(utr5, 'utr5', 1500)
    intron_df = load_df(intron, 'intron', 1500)
    cds_df = load_df(cds, 'cds', 1500)

    fig, ax = plt.subplots()
    fig, axs = plt.subplots(figsize=(12,12), ncols=2, nrows=3)

    data_all = [[(first_coding_exon_df['length'], 'First coding exon'),
                 (last_coding_exon_df['length'], 'Last coding exon')],
                [(utr3_df['length'], "3'UTR"), (utr5_df['length'], "5'UTR")],
                [(intron_df['length'], 'All Introns'), (cds_df['length'], 'CDS')]]

    for row in (0,1,2):
        for col in (0,1):
            data = data_all[row][col]
            print(row,col)
            sns.distplot(data[0], ax=axs[row, col],  kde=False, label=data[1], color='b')
            axs[row, col].set_title(data[1])#legend()
            axs[row, col].set_xlabel('')
    axs[0,0].set_ylabel('Frequency')
    axs[1,0].set_ylabel('Frequency')
    axs[2,0].set_ylabel('Frequency')
    fig.suptitle(title)
    plt.savefig('{}.png'.format(outprefix))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputdir', help='Input directory')
    parser.add_argument('title', help='Title')
    parser.add_argument('outprefix', help='Prefix of output file')
    args = parser.parse_args()
    plot_size_distribution(args.inputdir, args.title, args.outprefix)
