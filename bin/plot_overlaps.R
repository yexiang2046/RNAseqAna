#!/usr/bin/env Rscript

# Load required library
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    install.packages("VennDiagram")
}
library(VennDiagram)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript plot_overlaps.R <overlap_summary.csv>")
}

# Read the CSV file
csv_file <- args[1]
data <- read.csv(csv_file)

# Create output directory for plots
output_dir <- dirname(csv_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Process each row in the CSV
for (i in 1:nrow(data)) {
    row <- data[i,]
    
    # Create Venn diagram
    venn.plot <- draw.triple.venn(
        area1 = row$total1,
        area2 = row$total2,
        area3 = row$total3,
        n12 = row$overlap12,
        n13 = row$overlap13,
        n23 = row$overlap23,
        n123 = row$overlap123,
        category = c("File1", "File2", "File3"),
        fill = c("red", "blue", "green"),
        alpha = 0.5,
        label.col = "black",
        cex = 1,
        fontfamily = "sans",
        cat.col = c("red", "blue", "green"),
        cat.cex = 1,
        cat.fontfamily = "sans",
        print.mode = "raw"
    )
    
    # Save the plot
    output_file <- file.path(output_dir, paste0("venn_diagram_", i, ".pdf"))
    pdf(output_file, width = 8, height = 8)
    grid.draw(venn.plot)
    dev.off()
    
    cat(sprintf("Created Venn diagram: %s\n", output_file))
} 
