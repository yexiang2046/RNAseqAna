#!/usr/bin/env Rscript
# R script to analyze RAMPAGE data with CAGEr
# Usage: Rscript CAGEr_analysis.R <input_file> <output_file>
# input_file will contain bam file list and output_file will contain the output of the analysis as rds file
# input_file format: sample_name bam_file
# Load the required libraries

library(CAGEr)
library(GenomicRanges)
library(rtracklayer)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hkshv.encode.p13)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(argparse)

# Read the command line arguments
args <- commandArgs(trailingOnly = TRUE)


# Read the input file
input_file <- "iSLK_samples.txt"
output_label <- args[2]

input_data <- read.table(input_file, header = TRUE)
bam_files <- input_data[,2]

# Create CAGEexp object
cage_data <- CAGEexp(inputFiles = bam_files, inputFilesType = "bamPairedEnd", genomeName = "BSgenome.Hkshv.encode.p13", sampleLabels  = sub("/home/xiang/New_disk3/GEO_data/RAMPAGE/aligned/", "", sub( "Aligned.sortedByCoord.out.bam", "", basename(bam_files))))



# Quantify CAGE tags to identify TSS
iSLK <- getCTSS(cage_data, sequencingQualityThreshold = 10,
                mappingQualityThreshold = 20, 
                useMulticore = TRUE, nrCores = 12)


# Normalize tag counts
iSLK <- normalizeTagCount(iSLK, method = "simpleTpm") # You can choose other normalization methods

# Cluster CTSS to define consensus clusters (i.e., TSS regions)
iSLK <- clusterCTSS(iSLK, threshold = 1,
            method = "paraclu", useMulticore = TRUE,
            nrCores = 12) # The threshold can be adjusted



# Export results or perform additional analysis
tagClusters <- tagClustersGR(iSLK)


saveRDS(iSLK, file = "iSLK_CAGEr_RAMPAGE.rds")


# for BCBL1
# Read the input file
input_file <- "BCBL1_sample.txt"


input_data <- read.table(input_file, header = TRUE)
bam_files <- input_data[,2]

# Create CAGEexp object
cage_data <- CAGEexp(inputFiles = bam_files, inputFilesType = "bamPairedEnd", genomeName = "BSgenome.Hkshv.encode.p13", sampleLabels  = sub("/home/xiang/New_disk3/GEO_data/RAMPAGE/aligned/", "", sub( "Aligned.sortedByCoord.out.bam", "", basename(bam_files))))



# Quantify CAGE tags to identify TSS
BCBL1 <- getCTSS(cage_data, sequencingQualityThreshold = 10,
                mappingQualityThreshold = 20, 
                useMulticore = TRUE, nrCores = 12)



# Normalize tag counts
BCBL1 <- normalizeTagCount(BCBL1, method = "simpleTpm") # You can choose other normalization methods

# Cluster CTSS to define consensus clusters (i.e., TSS regions)
BCBL1 <- clusterCTSS(BCBL1, threshold = 1,
            method = "paraclu", useMulticore = TRUE,
            nrCores = 12) # The threshold can be adjusted



# Export results or perform additional analysis
tagClusters <- tagClustersGR(BCBL1)


saveRDS(BCBL1, file = "BCBL1_CAGEr_RAMPAGE.rds")

