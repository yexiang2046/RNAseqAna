# R script to count reads on each chromosome using bam file

# Load libraries
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)


# help message
if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  cat("Usage: Rscript Count_reads_on_each_chromosome.R bam_file1 bam_file2 ...\n")
  q()
}

# get bam file list from the command line
bam_files <- commandArgs(trailingOnly = TRUE)


# function to count reads on each chromosome
count_reads <- function(bam_file){
  # read bam file
  bam <- readGAlignmentPairs(bam_file)
  # get chromosome names
  chr_names <- seqnames(bam)
  # count reads on each chromosome
  chr_counts <- table(chr_names)
  # return chromosome counts
  return(chr_counts)
}

# apply function to each bam file
chr_counts <- lapply(bam_files, count_reads)

# convert list to data frame
chr_counts <- as.data.frame(do.call(rbind, chr_counts))

# save data frame to csv file
write.csv(chr_counts, "chr_counts.csv", row.names = TRUE)
