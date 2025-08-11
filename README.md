# RNA-seq Analysis Pipeline

A comprehensive Nextflow pipeline for RNA-seq data analysis, including quality control, alignment, feature counting, and differential expression analysis.

## Overview

This pipeline performs a complete RNA-seq analysis workflow:

1. **Read Trimming**: fastp for adapter trimming and quality filtering
2. **Genome Indexing**: STAR genome index generation
3. **Read Alignment**: STAR for RNA-seq read alignment
4. **Feature Counting**: featureCounts for gene expression quantification
5. **Differential Expression**: edgeR-based differential gene expression analysis
6. **Visualization**: PCA plots and DEG barplots

## Prerequisites

- **Nextflow** (>=21.04.0)
- **Docker** or container runtime
- **R** (>=4.0.0) with required packages:
  - edgeR
  - limma
  - pheatmap
  - RColorBrewer
  - optparse
  - tidyverse

## Pipeline Architecture

The pipeline is organized into modular processes:

- **TRIM**: Read trimming and quality filtering using fastp
- **STAR_INDEX**: Genome index generation for STAR
- **ALIGN**: RNA-seq read alignment using STAR
- **FEATURECOUNT**: Gene expression quantification
- **DE_ANALYSIS**: Differential expression analysis using edgeR

## Input Requirements

### Data folder with fastq.gz files
```
data/
├── sample1_R1_001.fastq.gz
├── sample1_R2_001.fastq.gz
├── sample2_R1_001.fastq.gz
├── sample2_R2_001.fastq.gz
└── ...
```

### Required Files
- **Genome FASTA**: Reference genome sequence file
- **GTF Annotation**: Gene annotation file (e.g., GENCODE)
- **Metadata**: Sample information file (tab-delimited)

### Metadata File Format
The `metadata.txt` file should be tab-delimited with:
- **SampleID**: Sample identifier
- **group**: Experimental group/condition

## Usage

### Basic Run
```bash
nextflow run rnaseq_pipeline.nf \
    --readspath "data" \
    --genome "path/to/genome.fa" \
    --gtf "path/to/annotation.gtf" \
    --metadata "metadata.txt" \
    --outdir "results" \
    --species human
```

### Parameter Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--readspath` | "data" | Directory containing FASTQ files |
| `--genome` | Required | Reference genome FASTA file |
| `--gtf` | Required | Gene annotation GTF file |
| `--metadata` | Required | Sample metadata file |
| `--outdir` | "results" | Output directory |
| `--cpus` | 12 | Number of CPU cores |
| `--ram` | 60000000000 | Memory limit (~60GB) |
| `--species` | "human" | Species for annotation |

### Output Structure
```
results/
├── trimmed/             # Trimmed FASTQ files
├── star_index/          # STAR genome index
├── aligned/             # Aligned BAM files
├── feature_counts/      # Gene expression counts
└── de_results/          # Differential expression results
    ├── *.csv           # DE analysis results
    ├── PCA_plot.png    # Principal component analysis
    └── DEG_barplot_*.png # Differential gene plots
```


