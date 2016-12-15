#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=2){
    stop('Insufficient parameters.\nUsage: ./create_regions_from_gencode.R <gff_file> <output_dir>')
}
suppressMessages(library(GenomicFeatures))
gencode_gff <- args[1]
output_dir <- args[2]

##TODO: remove redundant functions
create_df_names <- function(gr, filename, record.names=NULL){
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=record.names,
                   scores=c(rep('.', length(gr))),
                   strands=strand(gr))
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
}

create_df <- function(gr, filename){
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=c(rep('.', length(gr))),
                   scores=c(rep('.', length(gr))),
                   strands=strand(gr))
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
}

TxDb <- makeTxDbFromGFF(gencode_gff)


exons.data <- unique(unlist(exons(TxDb, columns='tx_name')))
introns.data <- unique(unlist(intronsByTranscript(TxDb)))
fiveUTRs.data <- unique(unlist(fiveUTRsByTranscript(TxDb)))
threeUTRs.data <- unique(unlist(threeUTRsByTranscript(TxDb)))
transcripts.data <- unique(unlist(transcripts(TxDb, columns=NULL)))
cds.data <- unique(unlist(cds(TxDb, columns=NULL)))
genes.data <- unique(unlist(genes(TxDb, columns='tx_name')))
promoters.data <- unique(unlist(promoters(TxDb, upstream=1000, downstream=1000)))

create_df(exons.data, file.path(output_dir, 'exons.bed'))
create_df(introns.data, file.path(output_dir, 'introns.bed'))
create_df(fiveUTRs.data, file.path(output_dir, '5UTRs.bed'))
create_df(threeUTRs.data, file.path(output_dir, '3UTRs.bed'))
create_df(transcripts.data, file.path(output_dir, 'transcripts.bed'))
create_df(cds.data, file.path(output_dir, 'cds.bed'))
create_df(promoters.data, file.path(output_dir, 'promoters.bed'))
create_df_names(genes.data, file.path(output_dir, 'genes.bed'), 
                names(mcols(genes.data)$tx_name))

