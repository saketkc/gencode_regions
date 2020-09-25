#!/usr/bin/env
"""
Extract lincRNA coordinates from GTF
"""
import sys
import GTF
import numpy as np
import pandas as pd


def main(GENCODE):
    gc = GTF.dataframe(GENCODE)
    gc.gene_id = gc.gene_id.replace(to_replace=r"\.[0-9]+", value="", regex=True)

    idx = (gc.feature == "transcript") & (gc.transcript_type == "lincRNA")
    lincRNA = gc.ix[idx, ["seqname", "start", "end", "gene_id", "gene_name", "strand"]]
    lincRNA.start = lincRNA.start.astype(int)
    lincRNA.end = lincRNA.end.astype(int)
    lincRNA.sort_values(by=["seqname", "start", "end"], inplace=True)
    lincRNA.to_csv("lincRNA.bed", sep="\t", header=False, index=False)

    idx = (gc.feature == "gene") & (gc.gene_type == "lincRNA")
    lincRNA = gc.ix[idx, ["seqname", "start", "end", "gene_id", "gene_name", "strand"]]
    lincRNA.start = lincRNA.start.astype(int)
    lincRNA.end = lincRNA.end.astype(int)
    lincRNA.sort_values(by=["seqname", "start", "end"], inplace=True)
    lincRNA.to_csv("lincRNA_genes.bed", sep="\t", header=False, index=False)


if __name__ == "__main__":
    main(sys.argv[1])
