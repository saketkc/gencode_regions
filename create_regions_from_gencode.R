library(GenomicFeatures)
args <- commandArgs(trailingOnly = TRUE)

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


exons <- unique(unlist(exons(TxDb, columns='tx_name')))
introns <- unique(unlist(intronsByTranscript(TxDb)))
fiveUTRs <- unique(unlist(fiveUTRsByTranscript(TxDb)))
threeUTRs <- unique(unlist(threeUTRsByTranscript(TxDb)))
transcripts <- unique(unlist(transcripts(TxDb, columns=NULL)))
CDS <- unique(unlist(cds(TxDb, columns=NULL)))
GENES <- unique(unlist(genes(TxDb, columns='tx_name')))

create_df(exons, file.path(output_dir, 'exons.bed'))
create_df(introns, file.path(output_dir, 'introns.bed'))
create_df(fiveUTRs, file.path(output_dir, '5UTRs.bed'))
create_df(threeUTRs, file.path(output_dir, '3UTRs.bed'))
create_df(transcripts, file.path(output_dir, 'transcripts.bed'))
create_df(CDS, file.path(output_dir, 'CDS.bed'))
create_df_names(GENES, file.path(output_dir, 'genes.bed'), 
                names(mcols(GENES)$tx_name))

