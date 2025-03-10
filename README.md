# RNASeq_tools
scripts for processing transcriptomic data from fastq files to DEGs and functional analysis

Functions

  1.  Trim with fastp and Align reads with STAR
  2.  generate bigwig file for Genome browsers
  3.  Call differentially expressed genes (DEGs)
  4.  Functional analyses on DEGs
      1. GSEA analyses using log2Foldchange ranked gene list
      2. ORA analyses with DEGs
  5.  Transcript quantification with Salmon


Dependencies

Command line tools
  1. samtools
  2. STAR
  3. salmon (>= 1.4.0)
  4. fastqc
  5. multiqc
 
Python tools
  1. deeptools
  
R tools
  1. tidyverse
  2. edgeR
  3. msigdbr
  4. fgsea
  5. clusterProfiler
  6. DESeq2
  7. tximport
  8. rmarkdown


Run
  1. trim 
    ```bash
    bash/fastp_trim.sh data _S1_L005_R1_001.fastq.gz _S1_L005_R2_001.fastq.gz
    ```
  2. align
    ```bash
    bash/star_align.sh /path/to//trimmed .R1.fastp.fastq.gz .R2.fastp.fastq.gz /path/to/star_idx /path/to/aligned /path/to//GRCh38.primary_assembly.genome.fa
    ```
  3. Salmon quantification (two options):

     a. Using Nextflow workflow:
     ```bash
     nextflow run workflows/salmon_quant.nf \
         --reads 'data/*_{1,2}.fastq.gz' \
         --transcriptome reference/transcriptome.fa \
         --genome reference/genome.fa \
         --metadata metadata.csv \
         --contrasts contrasts.csv \
         --outdir results
     ```

     b. Using bash script:
     ```bash
     bash workflows/run_salmon_pipeline.sh \
         -r data/ \
         -t reference/transcriptome.fa \
         -g reference/genome.fa \
         -o results \
         -m metadata.csv \
         -c contrasts.csv \
         -p 12
     ```

## Salmon Quantification Workflow

The Salmon workflow performs the following steps:

1. **Quality Control**:
   - FastQC analysis of raw reads
   - MultiQC report generation

2. **Salmon Indexing**:
   - Creates a decoy-aware index using both transcriptome and genome
   - Improves mapping accuracy by reducing false positives

3. **Quantification**:
   - Performs transcript quantification with bias correction
   - Includes GC and sequence bias correction
   - Automatic library type detection

4. **Differential Expression**:
   - EdgeR analysis of quantification results
   - Generates comprehensive report and visualizations

### Required Files

1. **Read Files**: Paired-end fastq files named `*_1.fastq.gz` and `*_2.fastq.gz`
2. **Reference Files**:
   - Transcriptome FASTA (e.g., GENCODE transcripts)
   - Genome FASTA (for decoy-aware index)
3. **Metadata**: CSV file with sample information
4. **Contrasts**: CSV file defining comparisons (optional)

# RNA-seq Analysis Pipeline

## Differential Expression Analysis with EdgeR

The pipeline includes a flexible differential expression analysis module using EdgeR that can handle both featureCounts and Salmon quantification outputs.

### Input Requirements

1. **Count Data**:
   - Either featureCounts output (tab-delimited count matrix)
   - Or Salmon quantification directory containing quant.sf files

2. **Metadata File** (CSV format):
   - Row names: sample IDs
   - Must include a 'condition' column specifying treatment groups
   ```csv
   sample_id,condition
   sample1,control
   sample2,control
   sample3,treatment
   sample4,treatment
   ```

3. **Contrasts File** (optional, CSV format):
   - Specifies which conditions to compare
   - Required columns: name, treatment, control
   ```csv
   name,treatment,control
   contrast1,treatment,control
   contrast2,treatment2,control
   ```

### Usage

The EdgeR analysis can be run in two ways:

1. **For featureCounts data**:
   ```bash
   Rscript bin/edger.r counts_matrix.txt metadata.csv contrasts.csv
   ```

2. **For Salmon quantification**:
   ```bash
   Rscript bin/edger.r salmon_quant_directory metadata.csv contrasts.csv
   ```

### Outputs

The analysis generates the following outputs in the `de_results` directory:

1. **CSV Results Files**:
   - One file per contrast
   - Contains: GeneID, logFC, logCPM, F-statistic, P-value, FDR
   - Named as: `<contrast_name>_edgeR_results.csv`

2. **Visualization**:
   - MA plots (`<contrast_name>_MA_plot.pdf`)
   - Volcano plots (`<contrast_name>_volcano_plot.pdf`)

3. **HTML Report** (`EdgeR_analysis_report.html`):
   - Analysis overview
   - Quality control plots (MDS, BCV)
   - Results summary for each contrast
   - Interactive tables and plots

### Quality Control

The analysis includes several QC steps:

1. **Filtering**: Removes low-count genes using EdgeR's `filterByExpr`
2. **Normalization**: Applies TMM normalization
3. **Dispersion Estimation**: Robust dispersion estimation
4. **Multiple Testing**: FDR correction using Benjamini-Hochberg method

### Key Features

- Automatic detection of input type (Salmon or featureCounts)
- Support for multiple contrasts
- Comprehensive QC visualizations
- Publication-ready plots
- Detailed HTML report generation
- Robust statistical analysis using EdgeR's GLM framework

### Dependencies

Required R packages:
- edgeR
- tximport (for Salmon data)
- tidyverse
- jsonlite
- rmarkdown
