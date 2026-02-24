#!/bin/bash

mkdir genome_files

cd genome_files

##get the fasta file 
wget https://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

##unzip fasta file
gzip -dk Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

#filter fasta file for chromosomes in vcf
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa \
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X \
   > Homo_sapiens.GRCh37.chr1-22_X.fa

##create indexed fasta file
samtools faidx Homo_sapiens.GRCh37.chr1-22_X.fa


##get gtf file
https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

##unzip gtf file
gzip -dk Homo_sapiens.GRCh37.87.gtf.gz


##filter gtf file for the same chromosomes as fasta file
awk '$1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X)$/ || /^#/' Homo_sapiens.GRCh37.87.gtf > filtered_chr1-22_X.gtf

cd ..
mkdir snp-files
cd snp-files

##get the vcf file
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz

##for the vcf file first edit the chr1 --> 1 as per the fasta file format

#Processing vcf files so as to keep the chr 1-22 and X information

##Extract the VCF header (This should not have the ##contigs)
bcftools view -h hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf > vcf_header.txt 

##Generate the contig header lines from the FASTA index (filtered fasta file) - this will have the ##contigs 1-22-X info
awk '{print "##contig=<ID="$1",length="$2">"}' ../genome_files/Homo_sapiens.GRCh37.chr1-22_X.fa.fai > ../genome_files/contigs_header.txt

##Insert the contig lines right after the ##fileformat line in the header: (creates a new file where the ##contigs info are added before the vcf contents)
awk '
  /^##fileformat/ {
    print $0
    system("cat ../genome_files/contigs_header.txt")
    next
  }
  {print $0}
' vcf_header.txt > new_header.txt

##Apply the new header to the VCF
bcftools reheader \
  -h new_header.txt \
  -o hg19_v0_1000G_phase1.snps.high_confidence.b37.withcontigs.vcf \
  hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf

##Verify if added
bcftools view -h hg19_v0_1000G_phase1.snps.high_confidence.b37.withcontigs.vcf | grep "^##contig"
Output should look like: 
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
...
##contig=<ID=X,length=155270560>

##Sanity check for chromosome filter names in fasta files and vcf files

##Check FASTA chromosomes

##Use the .fai index you already generated:
cut -f1 ../genome_files/Homo_sapiens.GRCh37.chr1-22_X.fa.fai

##Check VCF chromosomes
bcftools query -f '%CHROM\n' hg19_v0_1000G_phase1.snps.high_confidence.b37.withcontigs.vcf | sort -u

##Check for overlap between FASTA and VCF

##To make 100% sure:
comm -12 \
  <(cut -f1 ../genome_files/Homo_sapiens.GRCh37.chr1-22_X.fa.fai | sort) \
  <(bcftools query -f '%CHROM\n' hg19_v0_1000G_phase1.snps.high_confidence.b37.withcontigs.vcf | sort)

##After making sure the vcf file is formatted for the ##contigs and the chromosomes finalise the vcf file
##zip and index

bgzip -c hg19_v0_1000G_phase1.snps.high_confidence.b37.withcontigs.vcf \
  > hg19_starready.vcf.gz

tabix -p vcf hg19_starready.vcf.gz

##Quick check (should print 1..22 and X only):

bcftools query -f '%CHROM\n' hg19_starready.vcf.gz | sort -u



