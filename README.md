## Docker Images

Pre-built Docker images are available on Docker Hub:

| Environment | Tool | Pull Command |
|-------------|------|-------------|
| STAR | RNA-seq alignment | `docker pull ashviyer/hd-ase-star` |
| GATK | Variant calling | `docker pull ashviyer/hd-ase-gatk` |
| featureCounts | Read counting | `docker pull ashviyer/hd-ase-featurecount` |
| RSeQC | QC analysis | `docker pull ashviyer/hd-ase-rseqc` |
| pyNBS | Network clustering | `docker pull ashviyer/hd-ase-pynbs` |

### Run interactively:
docker run -it ashviyer/hd-ase-star

```mermaid
flowchart TD
    A[FASTA,GTF,VCF] --> B[Preprocessing files: processing_files.sh]
    B --> C[Create genome index files: STAR_create_genome_index.sh]
    C --> D[Fist STAR align: First_STAR_align.sh]
    D --> E[Prepare for GATK: run_gatk_steps.sh]
    E --> F[Generate vcf per sample: run_Haplotypecaller.sh]
    F --> G[Filter vcf for biallelic SNPs: filter_per_sample_vcf.sh]
    G --> H[Indev vcf file: index_per_sample_vcfs.sh]
    H --> I[Align trimmed reads with WASP correction: star_wasp_align.sh]
    I --> J[Process WASP corrected BAM files: add_rg_and_index_wasp.sh]
    O --> K[Get raw counts: run_featurecounts_wasp.sh]
    J --> L[ASECounter: ASEReadCounter.sh]
    M[Download fastq files: fetch_fastq_sra.sh] --> N[QC for reads: fastqc]
    N --> O[Trim adapter sequences: read_trimming_fastp.sh]
    O --> N
    O --> I[Align trimmed reads with WASP correction: star_wasp_align.sh]    
```

    ## Pipeline

| Step | Script | Environment | Input | Output |
|------|--------|-------------|-------|--------|
| 1. Alignment | `scripts/align.sh` | `star_env` | FASTQ | BAM |
| 2. QC | `scripts/qc.sh` | `rseqc_env` | BAM | QC report |
| 3. Read counting | `scripts/featurecount.sh` | `featureCount_env` | BAM | Count matrix |
| 4. Diff. expression | `scripts/deseq2.R` | `r_env` | Count matrix | DE genes |
| 5. Variant calling | `scripts/gatk.sh` | `gatk_env` | BAM | VCF |
| 6. Network clustering | `scripts/pynbs_parallel.py` | `pynbs_env` | Mutation matrix | Clusters |
| 7. Visualisation | `scripts/plots.R` | `r_env` | DE genes + Clusters | Figures |
