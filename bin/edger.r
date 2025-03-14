#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(edgeR)
    library(tximport)
    library(tidyverse)
    library(jsonlite)
})

# Function to check if input is from Salmon
is_salmon_input <- function(input_dir) {
    return(dir.exists(input_dir) && any(grepl("quant.sf$", list.files(input_dir, recursive=TRUE))))
}

# Function to read Salmon quantification data
read_salmon_data <- function(input_dir) {
    # Get all quant.sf files
    quant_files <- list.files(input_dir, pattern="quant.sf$", recursive=TRUE, 
                             full.names=TRUE)
    names(quant_files) <- basename(dirname(quant_files))
    
    # Read first quant file to get transcript info
    first_file <- read.table(quant_files[1], header=TRUE)
    tx_ids <- first_file$Name
    
    # Create tx2gene mapping (assuming transcript IDs contain gene IDs before .)
    # This handles common formats like ENST00000123.4 -> ENSG00000123
    tx2gene <- data.frame(
        TXNAME = tx_ids,
        GENEID = sub("\\.[0-9]+$", "", # Remove version numbers
                     sub("^([^.]+).*", "\\1", tx_ids)) # Get part before first dot
    )
    
    # Import with tximport using tx2gene mapping
    txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, txOut=FALSE)
    
    # Create DGEList object
    counts <- txi$counts
    y <- DGEList(counts=counts)
    
    # Add gene names
    y$genes <- data.frame(GeneID=rownames(counts))
    
    return(y)
}

# Function to read featureCounts data
read_featurecounts_data <- function(counts_file) {
    # Read counts data
    counts_data <- read.delim(counts_file, row.names=1)
    
    # Remove length and other non-count columns if present
    count_cols <- !grepl("Length|Chr|Start|End|Strand", colnames(counts_data))
    counts <- counts_data[, count_cols]
    
    # Create DGEList object
    y <- DGEList(counts=counts)
    y$genes <- data.frame(GeneID=rownames(counts))
    
    return(y)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: edger.r <input_path> <metadata_file> [contrasts_file]")
}

input_path <- args[1]
metadata_file <- args[2]
contrasts_file <- if(length(args) >= 3) args[3] else NULL

# Read metadata
metadata <- read.csv(metadata_file, row.names=1)

# Read count data based on input type
if (is_salmon_input(input_path)) {
    y <- read_salmon_data(input_path)
} else {
    y <- read_featurecounts_data(input_path)
}

# Ensure sample names match between counts and metadata
common_samples <- intersect(colnames(y), rownames(metadata))
y <- y[, common_samples]
metadata <- metadata[common_samples, ]

# Create output directory
dir.create("de_results", showWarnings=FALSE)

# Filter low count genes
keep <- filterByExpr(y, group=metadata$condition)
y <- y[keep, ]

# Calculate normalization factors
y <- calcNormFactors(y)

# Create design matrix
design <- model.matrix(~0 + condition, data=metadata)
colnames(design) <- levels(factor(metadata$condition))

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit model
fit <- glmQLFit(y, design)

# Function to perform DE analysis for a contrast
run_de_analysis <- function(contrast_name, condition1, condition2) {
    # Create contrast
    contrast <- makeContrasts(
        contrasts=paste0(condition1, "-", condition2),
        levels=design
    )
    
    # Perform test
    qlf <- glmQLFTest(fit, contrast=contrast)
    
    # Get results
    res <- topTags(qlf, n=Inf)
    res_df <- as.data.frame(res)
    
    # Add gene names if available
    if (!is.null(y$genes)) {
        res_df <- cbind(GeneID=rownames(res_df), res_df)
    }
    
    # Save results
    write.csv(res_df, 
              file=file.path("de_results", paste0(contrast_name, "_edgeR_results.csv")),
              row.names=FALSE)
    
    # Create MA plot
    pdf(file.path("de_results", paste0(contrast_name, "_MA_plot.pdf")))
    plotMD(qlf, main=paste("MA Plot -", contrast_name))
    abline(h=c(-1,1), col="blue")
    dev.off()
    
    # Create Volcano plot
    pdf(file.path("de_results", paste0(contrast_name, "_volcano_plot.pdf")))
    with(res_df, plot(logFC, -log10(PValue),
         pch=20, main=paste("Volcano Plot -", contrast_name)))
    with(subset(res_df, FDR < 0.05 & abs(logFC) > 1),
         points(logFC, -log10(PValue), pch=20, col="red"))
    dev.off()
    
    return(res_df)
}

# Process contrasts
if (!is.null(contrasts_file)) {
    contrasts <- read.csv(contrasts_file)
    for(i in 1:nrow(contrasts)) {
        run_de_analysis(
            contrasts$name[i],
            contrasts$treatment[i],
            contrasts$control[i]
        )
    }
} else {
    # Default comparison if no contrasts file provided
    conditions <- levels(factor(metadata$condition))
    if (length(conditions) >= 2) {
        run_de_analysis(
            "default_contrast",
            conditions[1],
            conditions[2]
        )
    }
}

# Generate summary report
writeLines('
---
title: "EdgeR Differential Expression Analysis Report"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Analysis Overview

- Input type: ', if(is_salmon_input(input_path)) "Salmon quantification" else "featureCounts",
'
- Number of samples analyzed: ', ncol(y), '
- Number of genes analyzed: ', nrow(y), '
- Conditions compared: ', paste(levels(factor(metadata$condition)), collapse=", "), '

## Quality Control

```{r qc_plots, echo=FALSE}
# MDS plot
plotMDS(y)

# BCV plot
plotBCV(y)
```

## Results Summary

```{r results_summary, echo=FALSE}
if (!is.null(contrasts_file)) {
    contrasts <- read.csv("', contrasts_file, '")
    for(i in 1:nrow(contrasts)) {
        results <- read.csv(file.path("de_results", 
                          paste0(contrasts$name[i], "_edgeR_results.csv")))
        
        cat("### ", contrasts$name[i], "\n\n")
        
        cat("- Total DE genes (FDR < 0.05):", sum(results$FDR < 0.05), "\n")
        cat("- Up-regulated:", sum(results$FDR < 0.05 & results$logFC > 0), "\n")
        cat("- Down-regulated:", sum(results$FDR < 0.05 & results$logFC < 0), "\n\n")
        
        cat("Top 10 differentially expressed genes:\n\n")
        print(knitr::kable(head(results, 10)))
        cat("\n\n")
        
        cat("MA Plot:\n\n")
        cat(sprintf("![MA Plot](%s)\n\n", 
            file.path("de_results", paste0(contrasts$name[i], "_MA_plot.pdf"))))
        
        cat("Volcano Plot:\n\n")
        cat(sprintf("![Volcano Plot](%s)\n\n", 
            file.path("de_results", paste0(contrasts$name[i], "_volcano_plot.pdf"))))
    }
}
```
', "report.Rmd")

# Render the report
rmarkdown::render("report.Rmd", output_file="EdgeR_analysis_report.html")





