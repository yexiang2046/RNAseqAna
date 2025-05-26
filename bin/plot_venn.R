#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    if (!require(VennDiagram)) install.packages("VennDiagram")
    if (!require(RColorBrewer)) install.packages("RColorBrewer")
    library(VennDiagram)
    library(RColorBrewer)
})

# Suppress VennDiagram info messages
futile.logger::flog.threshold(futile.logger::ERROR)

#' Create a publication-quality three-way Venn diagram
#' 
#' @param set1_size Size of set 1
#' @param set2_size Size of set 2
#' @param set3_size Size of set 3
#' @param n12 Number of overlaps between set 1 and 2
#' @param n23 Number of overlaps between set 2 and 3
#' @param n13 Number of overlaps between set 1 and 3
#' @param n123 Number of overlaps between all three sets
#' @param set_names Vector of three names for the sets
#' @param output_file Output PDF file path
#' @param colors Vector of three colors (optional)
#' @param title Main title for the diagram (optional)
plot_three_way_venn <- function(set1_size, set2_size, set3_size,
                               n12, n23, n13, n123,
                               set_names = c("Set 1", "Set 2", "Set 3"),
                               output_file = "venn_diagram.pdf",
                               colors = NULL,
                               title = "") {
    
    # Set default colors if not provided
    if (is.null(colors)) {
        colors <- brewer.pal(3, "Set2")
    }
    
    # Calculate semi-transparent fill colors
    transparent_colors <- sapply(colors, function(x) {
        rgb(t(col2rgb(x))/255, alpha = 0.5)
    })
    
    # Set up the PDF device
    pdf(output_file, width = 10, height = 8, useDingbats = FALSE)
    
    # Create a new plot
    plot.new()
    
    # Create the Venn diagram
    venn_plot <- draw.triple.venn(
        area1 = set1_size,
        area2 = set2_size,
        area3 = set3_size,
        n12 = n12,
        n23 = n23,
        n13 = n13,
        n123 = n123,
        category = set_names,
        col = colors,
        fill = transparent_colors,
        alpha = 0.50,
        cex = 1.2,
        cat.cex = 1.4,
        cat.col = colors,
        cat.dist = 0.05,
        cat.pos = c(-20, 20, 180),
        euler.d = TRUE,
        scaled = TRUE,
        lwd = 2,
        margin = 0.1
    )
    
    # Add title if provided
    if (title != "") {
        title(main = title, line = -1, cex.main = 1.5)
    }
    
    # Close the PDF device
    dev.off()
    
    # Print summary statistics
    cat("\nVenn Diagram Statistics:\n")
    cat("------------------------\n")
    cat(sprintf("%s: %d\n", set_names[1], set1_size))
    cat(sprintf("%s: %d\n", set_names[2], set2_size))
    cat(sprintf("%s: %d\n", set_names[3], set3_size))
    cat(sprintf("Overlap %s-%s: %d\n", set_names[1], set_names[2], n12))
    cat(sprintf("Overlap %s-%s: %d\n", set_names[2], set_names[3], n23))
    cat(sprintf("Overlap %s-%s: %d\n", set_names[1], set_names[3], n13))
    cat(sprintf("Common to all: %d\n", n123))
}

# Example usage (you can modify these numbers with your actual data)
if (!interactive()) {
    # Parse command line arguments if needed
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) == 0) {
        # Example usage with sample data
        plot_three_way_venn(
            area1 = 10576,
            area2 = 1014,
            area3 = 335,
            n12 = 633,
            n23 = 125,
            n13 = 179,
            n123 = 93,
            set_names = c("PKR_uninf", "PKR_wt", "PKR_de3l"),
            output_file = "PKR_ip_macs2_venn_diagram.pdf",
            colors = c("#4B89DC", "#37BC9B", "#F6BB42"),
            title = "Overlap between datasets"
        )
    } else {
        # You can add command-line argument parsing here if needed
        cat("Usage: Rscript plot_venn.R\n")
    }
} 


# plot J2 edger enrich venn diagram
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

# Read and process the data
files <- list.files(path = "J2ip_edger_deg", pattern = "*.csv", full.names = TRUE)
list_of_dfs <- lapply(files, read.csv)

# Process genes by logFC and padj
genes_by_logfc_padj <- list_of_dfs %>%
    map(~ .x %>%
        filter(logFC > 1 & padj < 0.05) %>%
        pull(ensembl))
names(genes_by_logfc_padj) <- c("J2_uninf", "J2_wt", "J2_de3l")

# Calculate set sizes and overlaps
set_sizes <- sapply(genes_by_logfc_padj, length)
overlaps <- list(
    n12 = length(intersect(genes_by_logfc_padj[[1]], genes_by_logfc_padj[[2]])),
    n23 = length(intersect(genes_by_logfc_padj[[2]], genes_by_logfc_padj[[3]])),
    n13 = length(intersect(genes_by_logfc_padj[[1]], genes_by_logfc_padj[[3]])),
    n123 = length(Reduce(intersect, genes_by_logfc_padj))
)

# Create the Venn diagram
plot_three_way_venn(
    set1_size = set_sizes[1],
    set2_size = set_sizes[2],
    set3_size = set_sizes[3],
    n12 = overlaps$n12,
    n23 = overlaps$n23,
    n13 = overlaps$n13,
    n123 = overlaps$n123,
    set_names = c("J2_uninf", "J2_wt", "J2_de3l"),
    output_file = "J2_ip_edger_enrich_venn_diagram.pdf",
    colors = c("#4B89DC", "#37BC9B", "#F6BB42"),
    title = "Overlap between J2 IP datasets"
)

