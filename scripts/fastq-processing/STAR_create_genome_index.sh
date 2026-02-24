#!/bin/bash

#intiate conda star environment
conda activate rnaseq_env

##Tip (prevents the old pthread issue): set a larger stack size in this shell before running STAR:

ulimit -s 65532

##RUN STAR 

##Create the hg19/ch37 genome index files 
STAR --runThreadN 30 \
  --runMode genomeGenerate \
  --genomeDir GrCh37_STAR/ \
  --genomeFastaFiles Homo_sapiens.GRCh37.chr1-22_X.fa \
  --sjdbGTFfile filtered_chr1-22_X.gtf \
  --varVCFfile snp-files/hg19_starready.vcf.gz
