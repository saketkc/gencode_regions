# gencode_regions
Extract 3'UTR, 5'UTR, CDS, Promoter, Genes from Gencode files

## Depedencies
   
- r>=3.2.1
- [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

## Run:

```{bash}
create_regions_from_gencode.R <path_to_GFF> <path_to_output_dir>
```

Will create `exons.bed, 3UTR.bed, 5UTR.bed, genes.bed, cds.bed` in `<output_dir>`


## Example

Download GFF(GRCh37, v25, comprehensive, CHR) from gencodegenes.org:
```wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz \
   && gunzip gencode.v25.annotation.gff3.gz```

Create regions:
```create_regions_from_gencode.R gencode.v25.annotation.gff3 /path/to/GRCh37/annotation```


