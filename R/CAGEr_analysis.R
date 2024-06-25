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
output_file <- args[2]

input_data <- read.table(input_file, header = TRUE)
bam_files <- input_data[,2]

# Create CAGEexp object
cage_data <- CAGEexp(inputFiles = bam_files, inputFilesType = "bamPairedEnd" genome = BSgenome.Hsapiens.UCSC.hg38, nthreads = 24)



