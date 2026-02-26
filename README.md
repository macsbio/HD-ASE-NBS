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
