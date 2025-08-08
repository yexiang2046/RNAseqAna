# RNA-seq Analysis Pipeline

A Nextflow pipeline for RNA-seq data analysis, including differential expression analysis using edgeR.

## Prerequisites

- Nextflow (>=21.04.0)
- R (>=4.0.0)
- Required R packages:
  - edgeR
  - limma
  - pheatmap
  - RColorBrewer
  - optparse
  - tidyverse

## Pipeline Overview

This pipeline performs:
1. Quality control and preprocessing of RNA-seq data
2. Differential expression analysis using edgeR
3. Visualization of results with heatmaps and plots

## Usage

### Step1: trim reads -> align -> count

```bash
nextflow run rnaseq_pipeline.nf \
    --cpus 8 \ # change base on your computer cpus
    --genome "genome.fasta" \
    --gtf "reference.gtf" \
    --outdir "results"
```

### Step2: differential gene expression analysis
```bash
docker run -w $(pwd) -v $(pwd):$(pwd) \
    xiang2019/rnaseq_renv:v1.0.1 \
    Rscript bin/edger.r \
    -c counts.txt \
    -m metadata.txt \
    -o output_dir \
    -g path/to/gtf \
    -s human
```

### Metadata File Format

The metadata.txt file should be tab-delimited with the following columns:
- SampleID
- group

