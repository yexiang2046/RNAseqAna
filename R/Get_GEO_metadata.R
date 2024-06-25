#!/usr/bin/env Rscript

# Ensure necessary packages are installed and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)

# Command-line interface function
fetch_geo_metadata <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("No GEO accession provided. Usage: Rscript this_script.R <GEO_accession>\n")
    quit(status = 1)
  }
  
  geo_accession <- args[1]
  
  # Load data from GEO
  cat("Fetching data for GEO accession:", geo_accession, "...\n")
  geo_data <- getGEO(geo_accession, destdir = tempdir())
  
  # Check if the data is a list (if multiple platforms are present)
  if (is.list(geo_data)) {
    geo_data <- geo_data[[1]]  # Using the first platform if multiple are present
  }
  

  # Optionally, return the full metadata if further analysis is required
  return(pData(phenoData(geo_data)))
}

# Execute the function if the script is called directly
if (interactive()) {
  cat("This script should be run from the command line.\n")
} else {
  write.table(fetch_geo_metadata(), file = paste0(geo_accession, "geo_metadata.txt"))
}
