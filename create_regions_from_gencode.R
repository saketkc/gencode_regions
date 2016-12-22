#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=2){
    stop('Insufficient parameters.\nUsage: ./create_regions_from_gencode.R <gff_file> <output_dir>')
}
suppressMessages(library(GenomicFeatures))
gencode_gff <- args[1]
output_dir <- args[2]
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
##TODO: remove redundant functions
create_df_names <- function(gr, filename, record.names=NULL){
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names= gsub('\\.[0-9]+', '', record.names),
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

exons.data <- unique(unlist(exons(TxDb, columns=c('GENEID', 'exon_rank'))))
introns.data <- unique(unlist(intronsByTranscript(TxDb)))
fiveUTRs.data <- unique(unlist(fiveUTRsByTranscript(TxDb, use.names=T)))
threeUTRs.data <- unique(unlist(threeUTRsByTranscript(TxDb, use.names=T)))
transcripts.data <- unique(unlist(transcripts(TxDb, columns=NULL)))
cds.data <- unique(unlist(cds(TxDb, columns=NULL)))
genes.data <- unique(unlist(genes(TxDb, columns='tx_name')))

## The exon_ranks can be ambiguous, we just take the consensus: mode of exon_ranks. This is not always correct, but then this is also not wrong.

create_df_names(exons.data, file.path(output_dir, 'exons.bed'), paste(mcols(exons.data)$GENEID, lapply(mcols(exons.data)$exon_rank, Mode), sep='__')  )
create_df(introns.data, file.path(output_dir, 'introns.bed'))
create_df_names(fiveUTRs.data, file.path(output_dir, '5UTRs.bed'), paste(mcols(fiveUTRs.data)$exon_name, mcols(fiveUTRs.data)$exon_rank, sep='__')  )
create_df_names(threeUTRs.data, file.path(output_dir, '3UTRs.bed'), paste(mcols(threeUTRs.data)$exon_name,mcols(threeUTRs.data)$exon_rank, sep='__')  )
create_df(transcripts.data, file.path(output_dir, 'transcripts.bed'))
create_df(cds.data, file.path(output_dir, 'cds.bed'))

## We still don't understand: What's a promoter?
promoters.length <- c(1000, 2000, 3000, 4000, 5000)
for (len in promoters.length){
    promoters.data <- unique(unlist(promoters(TxDb, upstream=len, downstream=len)))
    create_df(promoters.data, file.path(output_dir, paste('promoters', len, 'bed', sep='.')))
}
create_df_names(genes.data, file.path(output_dir, 'genes.bed'),
                names(mcols(genes.data)$tx_name))

