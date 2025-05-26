#!/usr/bin/env Rscript

# Function to install and load packages
install_and_load <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, 
                       repos = "https://cloud.r-project.org",
                       quiet = TRUE)
    }
    library(package, character.only = TRUE)
}

# Install and load required packages with error handling
tryCatch({
    # Install BiocManager if not already installed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    
    # List of required packages
    packages <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
    
    # Install and load each package
    invisible(sapply(packages, install_and_load))
    
}, error = function(e) {
    cat("Error installing/loading packages:", conditionMessage(e), "\n")
    quit(status = 1)
})

# Create the graph using DiagrammeR's DOT language directly
graph_dot <- "
digraph pipeline {
    rankdir=LR;
    node [style=filled, fontname=Arial];
    
    # Process nodes
    raw_bam [label='Raw BAM Files', shape=rectangle, fillcolor=lightyellow];
    trim [label='Trim', shape=rectangle, fillcolor=lightyellow];
    align [label='Align', shape=rectangle, fillcolor=lightyellow];
    macs2 [label='MACS2 Peak Calling', shape=rectangle, fillcolor=lightyellow];
    annotation [label='Peak Annotation', shape=rectangle, fillcolor=lightyellow];
    stats [label='Peak Statistics', shape=rectangle, fillcolor=lightyellow];
    viz [label='Visualization', shape=rectangle, fillcolor=lightyellow];
    
    # Output nodes
    peaks [label='NarrowPeak Files', shape=oval, fillcolor=lightgrey];
    anno_peaks [label='Annotated Peaks', shape=oval, fillcolor=lightgrey];
    plots [label='Distribution Plots', shape=oval, fillcolor=lightgrey];
    
    # Process edges
    raw_bam -> trim;
    trim -> align;
    align -> macs2;
    macs2 -> annotation;
    annotation -> stats;
    annotation -> viz;
    
    # Output edges [style=dashed]
    macs2 -> peaks [style=dashed];
    annotation -> anno_peaks [style=dashed];
    viz -> plots [style=dashed];
}
"

# Create a grViz object
grViz_obj <- DiagrammeR::grViz(graph_dot)

# Export as SVG
svg_content <- DiagrammeRsvg::export_svg(grViz_obj)
cat(svg_content, file = "pipeline_flowchart.svg")

# Convert SVG to PDF and PNG using rsvg
rsvg::rsvg_pdf("pipeline_flowchart.svg", "pipeline_flowchart.pdf", width = 12, height = 8)
rsvg::rsvg_png("pipeline_flowchart.svg", "pipeline_flowchart.png", width = 1200, height = 800)

# Print message
cat("\nPipeline flowchart has been generated:\n")
cat("- pipeline_flowchart.svg\n")
cat("- pipeline_flowchart.pdf\n")
cat("- pipeline_flowchart.png\n") 