#!/usr/bin/env
"""
Extract rRNA coordinates from GTF
"""
import sys
import GTF
import numpy as np
import pandas as pd


def main(GENCODE):
    gc = GTF.dataframe(GENCODE)
    gc.gene_id = gc.gene_id.replace(to_replace=r"\.[0-9]+", value="", regex=True)

    idx = (gc.feature == "transcript") & gc.transcript_type.str.contains("rRNA")
    rRNA = gc.ix[
        idx, ["seqname", "start", "end", "transcript_id", "gene_name", "strand"]
    ]
    rRNA.start = rRNA.start.astype(int)
    rRNA.end = rRNA.end.astype(int)
    rRNA.sort_values(by=["seqname", "start", "end"], inplace=True)
    rRNA.to_csv("rRNA_transcripts.bed", sep="\t", header=False, index=False)

    idx = (gc.feature == "gene") & gc.gene_type.str.contains("rRNA")
    rRNA = gc.ix[idx, ["seqname", "start", "end", "gene_id", "gene_name", "strand"]]
    rRNA.start = rRNA.start.astype(int)
    rRNA.end = rRNA.end.astype(int)
    rRNA.sort_values(by=["seqname", "start", "end"], inplace=True)
    rRNA.to_csv("rRNA_genes.bed", sep="\t", header=False, index=False)


if __name__ == "__main__":
    main(sys.argv[1])
