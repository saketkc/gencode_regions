#!/bin/bash
cat gencode.v25.annotation.gtf |awk '$0 !~ /gene_type "rRNA"/'|awk '$0 !~ /transcript_type "rRNA"/'|awk '$0 !~ /gene_type "Mt_rRNA"/'|awk '$0 !~ /transcript_type "Mt_rRNA"/' | awk '$0 !~ /gene_type "Mt_tRNA"/' | awk '$0 !~ /gene_type "tRNA"/' | awk '$0 !~ /transcript_type "tRNA"/' > gencode.v25.annotation.without_rRNA_tRNA.gtf

