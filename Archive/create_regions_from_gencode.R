#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=2){
    stop('Insufficient parameters.\nUsage: ./create_regions_from_gencode.R <gff_file> <output_dir>')
}
suppressMessages(library(GenomicFeatures))
suppressMessages(library(dplyr))
gencode_gff <- args[1]
output_dir <- args[2]

options(scipen=999)

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
  df <- arrange(df, seqnames, starts, ends, names)
  df <- unique(df)
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
  return(df)
}

create_df <- function(gr, filename){
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=c(rep('.', length(gr))),
                   scores=c(rep('.', length(gr))),
                   strands=strand(gr))
  df <- arrange(df, seqnames, starts, ends, names)
  df <- unique(df)
  write.table(df, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
  return(df)
}

unlist_obj <- function(obj) {
  unlist.obj <- unlist(obj, use.names=FALSE)
  names(unlist.obj) <- rep(names(obj), elementNROWS(obj))
  return (unlist.obj) 
}


TxDb <- makeTxDbFromGFF(gencode_gff)

transcripts.data <- transcripts(TxDb, columns=c("tx_name", "gene_id"))
anyDuplicated(elementMetadata(transcripts.data)$tx_name)

tx2gene <- unlist(elementMetadata(transcripts.data)$gene_id)
names(tx2gene) <- elementMetadata(transcripts.data)$tx_name

threeUTRs.data <- threeUTRsByTranscript(TxDb, use.names=T)
mcols(threeUTRs.data)$gene_id <- drop(transcripts.data$gene_id[match(names(threeUTRs.data),  transcripts.data$tx_name)])
names(threeUTRs.data) <- tx2gene[names(threeUTRs.data)]

fiveUTRs.data <- fiveUTRsByTranscript(TxDb, use.names=T)
mcols(fiveUTRs.data)$gene_id <- drop(transcripts.data$gene_id[match(names(fiveUTRs.data), transcripts.data$tx_name)])
names(fiveUTRs.data) <- tx2gene[names(fiveUTRs.data)]

cds.data <- cdsBy(TxDb, by='tx', use.names=T)
names(cds.data) <- tx2gene[names(cds.data)]

exons.data <- exonsBy(TxDb, by='tx', use.names=T)
names(exons.data) <- tx2gene[names(exons.data)]

introns.data <- intronsByTranscript(TxDb, use.names=T)
names(introns.data) <- tx2gene[names(introns.data)]

genes.data <- genes(TxDb, columns='tx_name')

all.introns <- unlist_obj(introns.data)
all.exons <- unlist_obj(exons.data)
all.cds <- unlist_obj(cds.data)
all.fiveUTRs <- unlist_obj(fiveUTRs.data)
all.threeUTRs <- unlist_obj(threeUTRs.data)

transcripts.data <- unlist(transcripts.data)
exons.data <- unlist(exons.data)
threeUTRs.data <- unlist(threeUTRs.data)
fiveUTRs.data <- unlist(fiveUTRs.data)
cds.data <- unlist(cds.data)
introns.data <- unlist(introns.data)
genes.data <- unlist(genes.data)


## The exon_ranks can be ambiguous, we just take the consensus: mode of exon_ranks. This is not always correct, but then this is also not wrong.
exons.df <- create_df_names(exons.data, file.path(output_dir, 'exons_rank.bed'), paste(names(exons.data), lapply(mcols(exons.data)$exon_rank, Mode), sep='__')  )
fiveutrs.df <- create_df_names(fiveUTRs.data, file.path(output_dir, '5UTRs_rank.bed'), paste(names(fiveUTRs.data), mcols(fiveUTRs.data)$exon_rank, sep='__'))
threeutrs.df <- create_df_names(threeUTRs.data, file.path(output_dir, '3UTRs_rank.bed'), paste(names(threeUTRs.data), mcols(threeUTRs.data)$exon_rank, sep='__') )

all_exons.df <- create_df_names(all.exons, file.path(output_dir, 'exons.bed'), names(all.exons))
all_introns.df <- create_df_names(all.introns, file.path(output_dir, 'introns.bed'), names(all.introns))
cds.df <- create_df_names(all.cds, file.path(output_dir, 'cds.bed'), names(all.cds))
cds.longestORF.df <- group_by(cds.df, seqnames, starts, names, scores, strands)
cds.codons <- summarise(cds.longestORF.df, endsM = max(ends))

write.table(cds.codons[c('seqnames', 'starts', 'endsM', 'names', 'scores', 'strands')], 
            file=file.path(output_dir, 'cds_maxORF.bed'), 
            quote=F, sep='\t', row.names=F, col.names=F)

all_fiveutrs.df <- create_df_names(all.fiveUTRs, file.path(output_dir, '5UTRs.bed'), names(all.fiveUTRs))
all_threeutrs.df <- create_df_names(all.threeUTRs, file.path(output_dir, '3UTRs.bed'), names(all.threeUTRs))
all_genes.df <- create_df_names(genes.data, file.path(output_dir, 'genes.bed'), names(mcols(genes.data)$tx_name))
all_transcripts.df <- create_df(transcripts.data, file.path(output_dir, 'transcripts.bed'))

tx2gene.df <- as.data.frame(tx2gene)
rownames(tx2gene.df) <- gsub('\\.[0-9]+', '', rownames(tx2gene.df))
tx2gene.df$tx <- rownames(tx2gene.df)
names(tx2gene.df) <- sub("tx2gene", "gene", names(tx2gene.df))
tx2gene.df$gene <- gsub('\\.[0-9]+', '', tx2gene.df$gene)
tx2gene <- tx2gene[c('tx', 'gene')]
write.table(tx2gene.df, file=file.path(output_dir, 'tx2gene.bed'), quote=F, sep='\t', row.names=F, col.names=F)

## We still don't understand: What's a promoter?
promoters.length <- c(1000, 2000, 3000, 4000, 5000)
for (len in promoters.length){
    promoters.data <- unlist(promoters(TxDb, upstream=len, downstream=len))
    promoters.df <- create_df(promoters.data, file.path(output_dir, paste('promoters', len, 'bed', sep='.')))
}

