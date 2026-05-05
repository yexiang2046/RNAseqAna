# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RNAseqAna is an RNA-seq analysis pipeline built with Nextflow. It runs each stage in a Docker container and supports local and AWS Batch execution. The pipeline produces a gene-level count matrix (featureCounts) from raw FASTQ reads.

## Common Commands

### Setup (Red Hat / CentOS / Amazon Linux 2023)

```bash
bash prepare.sh                        # installs Java 17, Nextflow, Docker
bash install_for_ami2023linux_aws.sh   # full setup including Docker CE on AWS AMI 2023
```

### Running the Pipeline

```bash
# Paired-end (default) — place FASTQs in ./data/ matching *{1,2}*.fastq.gz
nextflow run main.nf

# Single-end
nextflow run main.nf --single_end true

# Override output directory
nextflow run main.nf --outdir my_results/

# AWS Batch
nextflow run main.nf -profile awsbatch

# Resume a failed run
nextflow run main.nf -resume
```

### Running Tests

```bash
cd tests && bash run_tests.sh
# or directly:
nextflow run tests/test_pipeline.nf
```

## Pipeline Architecture

`main.nf` defines the `RNASEQ` workflow. Stages run in order:

| Stage | Module | Tool | Container |
|-------|--------|------|-----------|
| (optional) FASTQC | `modules/fastqc.nf` | FastQC | `xiang2019/rnaseq_cmd:v1.0.0` |
| TRIM | `modules/fastp_trim.nf` | fastp | `staphb/fastp:0.24.0` |
| STAR_INDEX | `modules/star_align.nf` | STAR | `alexdobin/star:2.6.1d` |
| ALIGN | `modules/star_align.nf` | STAR | `alexdobin/star:2.6.1d` |
| FEATURECOUNT | `modules/featurecount.nf` | featureCounts | `xiang2019/rnaseq_cmd:v1.0.0` |
| MULTIQC | `modules/multiqc.nf` | MultiQC | `multiqc/multiqc:pdf-v1.34` |

FASTQC is imported but commented out in the workflow. Enable it by uncommenting the `FASTQC(read_pairs_ch)` call in `main.nf`.

ALIGN uses `maxForks 1` — only one STAR alignment runs at a time to avoid memory contention.

MULTIQC runs last and aggregates three QC sources: fastp JSON reports (from TRIM), STAR `Log.final.out` files (from ALIGN), and the featureCounts summary. Output lands in `<outdir>/multiqc_report/`.

### Single-end vs Paired-end

`params.single_end` (default: `false`) controls read mode across all stages:

- **Input channel** (`main.nf`): `fromFilePairs` glob `*{1,2}*.fastq.gz` for paired-end; `fromPath` glob `*.fastq.gz` for single-end.
- **TRIM** (`fastp_trim.nf`): paired-end adds `--in2`, `--out2`, `--detect_adapter_for_pe`, `--correction`; single-end omits them.
- **ALIGN** (`star_align.nf`): paired-end passes both files to `--readFilesIn`; single-end passes only `reads[0]`.
- **FEATURECOUNT** (`featurecount.nf`): paired-end adds `-p --countReadPairs -B`; single-end omits them.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--gtf` | `null` | Path to annotation GTF (falls back to `params.projectDir` default) |
| `--single_end` | `false` | Set `true` for single-end reads |
| `params.cpus` | `12` | CPUs passed to STAR_INDEX |
| `params.ram` | `60 GB` | RAM hint for STAR |

Resource allocations in `nextflow.config`:
- STAR_INDEX / ALIGN: 16 CPUs, 64 GB
- FEATURECOUNT: 8 CPUs (process-level default)
- Default for all other processes: 8 CPUs, 32 GB

### Execution Profiles (`nextflow.config`)

- `standard` — local Docker, 8 CPUs / 64 GB
- `awsbatch` — AWS Batch executor with Fusion filesystem (S3 access via `fusion.enabled = true`)

## fastp Trimming Settings

Applied to all reads regardless of mode: Q10 quality threshold, 40% unqualified-base limit, max 5 Ns, 36 bp minimum length, poly-G/X trimming, 3′ sliding-window quality trimming (window 4, mean Q20).

## Test Infrastructure

`tests/test_pipeline.nf` generates synthetic BAM, GTF, FASTA, and RMSK files, then runs pipeline processes against them. `tests/run_tests.sh` drives the run and checks that expected output files exist.
