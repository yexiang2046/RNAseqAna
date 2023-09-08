# RNASeq_tools
scripts for processing transcriptomic data from fastq files to DEGs and functional analysis

Funtions

  1.  Align reads with STAR
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
