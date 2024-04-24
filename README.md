# RNASeq_tools
scripts for processing transcriptomic data from fastq files to DEGs and functional analysis

Funtions

  1.  Trim with fastp and Align reads with STAR
  2.  generate bigwig file for Genome browsers
  3.  Call differentially expressed genes (DEGs)
  4.  Functional analyses on DEGs
      1. GSEA analyses using log2Foldchange ranked gene list
      2. ORA analyses with DEGs


Dependencies

Command line tools
  1. samtools
  2. STAR
 
Python tools
  1. deeptools
  
R tools
  1. tidyverse
  2. edgeR
  3. msigdbr
  4. fgsea
  5. clusterProfiler
  6. DESeq2


Run
  1. trim 
	```shell
		bash/fastp_trim.sh data _S1_L005_R1_001.fastq.gz _S1_L005_R2_001.fastq.gz
	```
  2. align
	```shell
		bash/star_align.sh /path/to//trimmed .R1.fastp.fastq.gz .R2.fastp.fastq.gz /path/to/star_idx /path/to/aligned /path/to//GRCh38.primary_assembly.genome.fa
	```
