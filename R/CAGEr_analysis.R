#!/usr/bin/env Rscript
# R script to analyze RAMPAGE data with CAGEr
# Usage: Rscript CAGEr_analysis.R <input_file> <output_file>
# input_file will contain bam file list and output_file will contain the output of the analysis as rds file
# input_file format: sample_name bam_file
# Load the required libraries

library(CAGEr)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(argparse)

# Read the command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Read the input file
input_file <- args[1]
output_label <- args[2]

input_data <- read.table(input_file, header = TRUE)
bam_files <- input_data[,2]

# Create CAGEexp object
cage_data <- CAGEexp(inputFiles = bam_files, inputFilesType = "bamPairedEnd" genome = BSgenome.Hsapiens.UCSC.hg38, nthreads = 24)

# Set sequencing protocol to RAMPAGE
seqProt(cage_data) <- "RAMPAGE"

# Specify the genome if not set during initialization
genome(cage_data) <- "BSgenome.Hsapiens.UCSC.hg38.KSHV" # adjust as needed

# Quantify CAGE tags to identify TSS
getCTSS(cage_data)
quantifyCTSS(cage_data)

# Normalize tag counts
normalizeTagCount(cage_data, method = "quantile") # You can choose other normalization methods

# Cluster CTSS to define consensus clusters (i.e., TSS regions)
clusterCTSS(cage_data, threshold = 100) # The threshold can be adjusted

# Aggregate TSSs into tag clusters
aggregateTagClusters(cage_data)

# Export results or perform additional analysis
tagClusters <- getTagClusters(cage_data)

# Optionally, output the tag clusters to a file
write.table(tagClusters, file = paste0(output_label,"_tagClusters.txt"), sep = "\t", quote = FALSE)
saveRDS(cage_data, file = paste0(output_label,"_cage_data.rds")


