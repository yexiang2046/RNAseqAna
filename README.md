# RNA-seq Analysis Pipeline

A comprehensive Nextflow pipeline for RNA-seq data analysis, including quality control, alignment, and feature counting. Downstream differential expression and functional enrichment analysis are run independently using the provided R scripts.

## Overview

### Nextflow Pipeline (`main.nf`)

1. **Read Trimming**: fastp for adapter trimming and quality filtering
2. **Genome Indexing**: STAR genome index generation
3. **Read Alignment**: STAR for RNA-seq read alignment
4. **Feature Counting**: featureCounts for gene expression quantification
5. **Quality Control**: MultiQC aggregated report

### Downstream Analysis (independent R scripts)

6. **Differential Expression**: edgeR-based differential gene expression (`bin/edger.r`)
7. **Functional Analysis**: GO, KEGG, Reactome, and GSEA enrichment analysis (`bin/functional_analysis.r`)
8. **Heatmap**: DEG expression heatmap across samples (`bin/heatmap.r`)

## Prerequisites

- **Nextflow** (>=21.04.0)
- **Docker** or container runtime
- **Docker** image `xiang2019/rnaseq_renv:v1.0.2` for all downstream R scripts (edgeR, functional analysis, heatmap)

## Pipeline Architecture

The pipeline is organized into modular processes:

| Stage | Module | Tool | Container |
|-------|--------|------|-----------|
| TRIM | `modules/fastp_trim.nf` | fastp | `staphb/fastp:0.24.0` |
| STAR_INDEX | `modules/star_align.nf` | STAR | `quay.io/biocontainers/star:2.7.11b` |
| ALIGN | `modules/star_align.nf` | STAR | `quay.io/biocontainers/star:2.7.11b` |
| FEATURECOUNT | `modules/featurecount.nf` | featureCounts | `xiang2019/rnaseq_cmd:v1.0.0` |
| MULTIQC | `modules/multiqc.nf` | MultiQC | `multiqc/multiqc:pdf-v1.34` |

FASTQC is imported but disabled by default. Enable it by uncommenting `FASTQC(read_pairs_ch)` in `main.nf`.

ALIGN uses `maxForks 1` — only one STAR alignment runs at a time to avoid memory contention.

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
- **Genome FASTA**: Reference genome sequence file (`*.genome.fa` in project directory)
- **GTF Annotation**: Gene annotation file (e.g., GENCODE)
- **Metadata**: Tab-delimited sample information file for downstream DE and functional analysis

### Metadata File Format
The `metadata.txt` file should be tab-delimited with:
- **SampleID**: Sample identifier matching BAM file names
- **group**: Experimental group/condition

## Usage

### Running the Nextflow Pipeline

```bash
# Paired-end (default) — FASTQs in ./data/ matching *{1,2}*.fastq.gz
nextflow run main.nf

# Single-end
nextflow run main.nf --single_end true

# Use pre-built STAR index
nextflow run main.nf --star_index /path/to/star_index

# AWS Batch
nextflow run main.nf -profile awsbatch

# Resume a failed run
nextflow run main.nf -resume
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--gtf` | project default | Gene annotation GTF file |
| `--single_end` | `false` | Set `true` for single-end reads |
| `--star_index` | `null` | Path to pre-built STAR index (skips STAR_INDEX if set) |
| `--cpus` | `12` | CPUs for STAR indexing |
| `--ram` | `60 GB` | RAM hint for STAR |

### Differential Expression Analysis

Run after the Nextflow pipeline produces a count matrix:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/edger.r \
        -c /data/results/feature_counts/counts.txt \
        -m /data/metadata.txt \
        -g /data/annotation.gtf \
        -o /data/de_results \
        -s human
```

| Option | Description |
|--------|-------------|
| `-c` | featureCounts output file |
| `-m` | Metadata file (tab-delimited, columns: SampleID, group) |
| `-g` | GTF annotation file |
| `-o` | Output directory |
| `-s` | Species: `human` or `mouse` |
| `-p` | (Optional) Metadata column for paired/blocking factor |
| `-f` | (Optional) Comma-separated metadata columns to combine with `group` |

### Functional Enrichment Analysis

Run on each DEG file produced by `edger.r`:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/functional_analysis.r \
        /data/de_results/DEG_<comparison>.csv \
        /data/functional_analysis_<comparison> \
        human \
        <comparison_name>
```

Pass `FALSE` as a 5th argument to use FDR < 0.05 only (default requires FDR < 0.05 **and** |logFC| > 1):

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/functional_analysis.r \
        /data/de_results/DEG_KO_vs_WT.csv \
        /data/functional_analysis_KO_vs_WT \
        human KO_vs_WT FALSE
```

### Heatmap

Run on DEG files and normalized counts from `edger.r`:

```bash
# Single comparison
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/heatmap.r \
        -d /data/de_results/DEG_KO_vs_WT.csv \
        -n /data/de_results/normalized_counts_CPM.csv \
        -m /data/metadata.txt \
        -o /data/heatmaps

# Union of DEGs from multiple comparisons
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/heatmap.r \
        -d /data/de_results/DEG_KO_vs_WT.csv,/data/de_results/DEG_OE_vs_WT.csv \
        -n /data/de_results/normalized_counts_CPM.csv \
        -m /data/metadata.txt \
        --prefix KO_OE_union \
        -o /data/heatmaps

# Multiple counts files (e.g. separate experiments), with sample ordering and gaps
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/heatmap.r \
        -d /data/de_results/DEG_KO_vs_WT.csv \
        -n /data/exp1/normalized_counts_CPM.csv,/data/exp2/normalized_counts_CPM.csv \
        --samples "WT_1,WT_2,WT_3,KO_1,KO_2,KO_3" \
        --gap_samples "WT_3" \
        -o /data/heatmaps

# Stricter thresholds, FDR only (no logFC filter), cap at 50 genes
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(pwd):/data" \
    xiang2019/rnaseq_renv:v1.0.2 \
    Rscript /data/bin/heatmap.r \
        -d /data/de_results/DEG_KO_vs_WT.csv \
        -n /data/de_results/normalized_counts_CPM.csv \
        --fdr 0.01 --use_logfc FALSE --max_genes 50 \
        -o /data/heatmaps
```

| Option | Default | Description |
|--------|---------|-------------|
| `-d` / `--deg` | required | Comma-separated DEG CSV file(s); significant genes are unioned |
| `-n` / `--norm_counts` | required | Comma-separated normalized counts CPM CSV file(s); merged by `gene_id` |
| `-m` / `--metadata` | — | Metadata file for column annotation (SampleId + group columns) |
| `-o` / `--output` | `.` | Output directory |
| `-f` / `--fdr` | `0.05` | FDR threshold |
| `-l` / `--logfc` | `1.0` | Minimum \|logFC\| threshold |
| `--use_logfc` | `TRUE` | `FALSE` to filter by FDR only |
| `--max_genes` | `100` | Max genes to plot, ranked by max \|logFC\| across comparisons |
| `--samples` | all | Comma-separated sample IDs to include, in display order |
| `--gap_samples` | — | Comma-separated sample IDs after which to draw a column gap; disables column clustering |
| `--prefix` | from filename | Output file prefix |

## Functional Analysis Details

`bin/functional_analysis.r` runs the following analyses for each comparison:

| Analysis | Method | Database |
|----------|--------|----------|
| GO Biological Process | Over-representation (ORA) | Gene Ontology |
| GO Molecular Function | Over-representation (ORA) | Gene Ontology |
| GO Cellular Component | Over-representation (ORA) | Gene Ontology |
| KEGG Pathway | Over-representation (ORA) | KEGG |
| Reactome Pathway | Over-representation (ORA) | Reactome |
| GSEA | `fgseaMultilevel` ranked by logFC | MSigDB HALLMARK |
| ORA HALLMARK | Over-representation (ORA) | MSigDB HALLMARK |

Highly overlapping gene sets (Jaccard similarity > 70%) are filtered. Both full (`_full.csv`) and filtered results are saved for each analysis.

Supported species: `human` (hg38, `org.Hs.eg.db`) and `mouse` (mm10, `org.Mm.eg.db`).

## Output Structure

```
results/
├── trimmed/             # Trimmed FASTQ files and fastp JSON reports
├── star_index/          # STAR genome index
├── aligned/             # Aligned BAM files and STAR logs
├── feature_counts/      # Gene expression count matrix (counts.txt)
└── multiqc_report/      # Aggregated QC (fastp + STAR + featureCounts)

de_results/              # From edger.r
├── *.csv                # Full DE results per comparison
├── DEG_*.csv            # Significant DEGs per comparison
└── PCA_plot.png         # Principal component analysis

functional_analysis_<comparison>/    # From functional_analysis.r
├── GO_BP_<comparison>.csv           # GO BP (redundancy-filtered)
├── GO_BP_<comparison>_full.csv      # GO BP (full)
├── GO_BP_<comparison>.pdf
├── GO_MF_<comparison>.csv / .pdf
├── GO_CC_<comparison>.csv / .pdf
├── KEGG_<comparison>.csv / .pdf
├── Reactome_<comparison>.csv / .pdf
├── GSEA_HALLMARK_<comparison>.csv / .pdf
└── ORA_HALLMARK_<comparison>.csv / .pdf

heatmaps/                            # From heatmap.r
├── heatmap_<prefix>.pdf
└── heatmap_<prefix>.png
```

## Setup

```bash
bash prepare.sh                        # installs Java 17, Nextflow, Docker
bash install_for_ami2023linux_aws.sh   # full setup on AWS AMI 2023
```

## Running Tests

```bash
cd tests && bash run_tests.sh
# or directly:
nextflow run tests/test_pipeline.nf
```
