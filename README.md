# gencode_regions

Extract 3'UTR, 5'UTR, CDS, Promoter, Genes from GTF files.

# Using Python

## Dependencies
   
- Python>=3
- [gffutils](https://github.com/daler/gffutils)
- [pybedtools](https://daler.github.io/pybedtools)

## Notebooks

- [BDGP6](notebooks/BDGP6.ipynb)
- [GRCg6](notebooks/GRCg6.ipynb)
- [GRCz10](notebooks/GRCz10.ipynb)
- [GRCz11](notebooks/GRCz11.ipynb)
- [GRCh38](notebooks/hg38-v96.ipynb)
- [MG1655](notebooks/MG1655.ipynb)
- [GRCm10](notebooks/mm10.ipynb)
- [Mmul8](notebooks/Mmul8.ipynb)
- [panTro3](notebooks/panTro3.ipynb)
- [Rnor6.0](notebooks/Rnor6.0.ipynb)
- [sacCerR64](notebooks/sacCerR64.ipynb)
- [WBcel235](notebooks/WBcel235.ipynb)

# Using R

## Dependencies
   
- r>=3.2.1
- [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

## Run

```{bash}
./create_regions_from_gencode.R <path_to_GFF/GTF> <path_to_output_dir>
```

Will create `exons.bed, 3UTR.bed, 5UTR.bed, genes.bed, cds.bed` in `<output_dir>`


## Example

- Download GFF/GTF(GRCh37, v25, comprehensive, CHR) from gencodegenes.org:

```{bash}
   wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz \
   && gunzip gencode.v25.annotation.gff3.gz
```

- Create regions:

```{bash}
./create_regions_from_gencode.R gencode.v25.annotation.gff3 /path/to/GRCh37/annotation
```




## First exons, Last exons

We use [`GenePred`](https://genome.ucsc.edu/FAQ/FAQformat#format9) format to make the process a bit simple.

 - Download [gtfToGenePred](http://hgdownload.cse.ucsc.edu/admin/exe/)
 - Convert gtf to GenePred:
 
     ```{bash}
     gtfToGenePred gencode.v25.annotation.gtf gencode.v25.annotation.genepred
     ```
     
 - Extract `first exons`:
 
     ```{bash}
     python genepred_to_bed.py --first_exon gencode.v25.annotation.genepred
     ```
     
 - Extract `last exons`:
 
     ```{bash}
     python genepred_to_bed.py --last_exon gencode.v25.annotation.genepred
     ```





## Confused about exons and UTRs?

This should be helpful:
![img](images/eukaryotic_regulation.png)

[Source: Wikipedia](https://en.wikipedia.org/wiki/File:Gene_structure_eukaryote_2_annotated.svg)

or probably this:

![img](images/transcription_elements.jpg)

[Source: Biostar](https://www.biostars.org/p/47022/)



