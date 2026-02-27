```mermaid
flowchart TD
  %% Inputs
  REF[FASTA / GTF / (known SNP VCF)] --> PREP[Preprocess reference files<br/>processing_files.sh]
  PREP --> IDX[Build STAR genome index<br/>STAR_create_genome_index.sh]

  FASTQ[Download FASTQ<br/>fetch_fastq_sra.sh] --> QC1[QC (FastQC)]
  QC1 --> TRIM[Trim adapters<br/>read_trimming_fastp.sh]
  TRIM --> QC2[QC after trimming (FastQC) (optional)]

  %% Standard alignment + gene counts (expression)
  IDX --> STAR_STD[STAR align (standard)<br/>First_STAR_align.sh]
  TRIM --> STAR_STD
  QC2 --> STAR_STD
  STAR_STD --> FC[Gene counts (featureCounts)<br/>run_featurecounts.sh]

  %% Variant calling to get per-sample VCF for WASP
  STAR_STD --> GATK[Prepare BAM for GATK<br/>run_gatk_steps.sh]
  GATK --> HC[HaplotypeCaller per sample<br/>run_Haplotypecaller.sh]
  HC --> VCF_FILT[Filter to biallelic SNPs<br/>filter_per_sample_vcf.sh]
  VCF_FILT --> VCF_IDX[Index per-sample VCF<br/>index_per_sample_vcfs.sh]

  %% WASP alignment + ASE
  IDX --> STAR_WASP[STAR align with WASP correction<br/>star_wasp_align.sh]
  TRIM --> STAR_WASP
  VCF_IDX --> STAR_WASP
  STAR_WASP --> BAM_PROC[Process WASP BAM (add RG + index)<br/>add_rg_and_index_wasp.sh]
```
  BAM_PROC --> ASE[ASEReadCounter<br/>ASEReadCounter.sh]

