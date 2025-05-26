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
digraph rnaseq_pipeline {
    rankdir=LR;
    node [style=filled, fontname=Arial];
    
    # Input nodes
    fastq [label='Raw FASTQ Files', shape=oval, fillcolor=lightgrey];
    gtf [label='GTF Annotation', shape=oval, fillcolor=lightgrey];
    
    # Process nodes
    trim [label='Quality Trimming', shape=rectangle, fillcolor=lightyellow];
    align [label='Alignment (STAR)', shape=rectangle, fillcolor=lightyellow];
    sort [label='BAM Sorting', shape=rectangle, fillcolor=lightyellow];
    index [label='BAM Indexing', shape=rectangle, fillcolor=lightyellow];
    count [label='Feature Counting', shape=rectangle, fillcolor=lightyellow];
    diff [label='Differential Expression', shape=rectangle, fillcolor=lightyellow];
    enrich [label='Pathway Enrichment', shape=rectangle, fillcolor=lightyellow];
    
    # Output nodes
    bam [label='Sorted BAM Files', shape=oval, fillcolor=lightgrey];
    counts [label='Count Matrix', shape=oval, fillcolor=lightgrey];
    deg [label='DEG Results', shape=oval, fillcolor=lightgrey];
    pathways [label='Enrichment Results', shape=oval, fillcolor=lightgrey];
    
    # Process edges
    fastq -> trim;
    trim -> align;
    align -> sort;
    sort -> index;
    sort -> bam;
    bam -> count;
    gtf -> count;
    count -> counts;
    counts -> diff;
    diff -> deg;
    deg -> enrich;
    enrich -> pathways;
    
    # Output edges [style=dashed]
    bam -> bam [style=dashed];
    counts -> counts [style=dashed];
    deg -> deg [style=dashed];
    pathways -> pathways [style=dashed];
}
"

# Create a grViz object
grViz_obj <- DiagrammeR::grViz(graph_dot)

# Export as SVG
svg_content <- DiagrammeRsvg::export_svg(grViz_obj)
cat(svg_content, file = "rnaseq_pipeline_flowchart.svg")

# Convert SVG to PDF and PNG using rsvg
rsvg::rsvg_pdf("rnaseq_pipeline_flowchart.svg", "rnaseq_pipeline_flowchart.pdf", width = 12, height = 8)
rsvg::rsvg_png("rnaseq_pipeline_flowchart.svg", "rnaseq_pipeline_flowchart.png", width = 1200, height = 800)

# Print message
cat("\nRNA-seq pipeline flowchart has been generated:\n")
cat("- rnaseq_pipeline_flowchart.svg\n")
cat("- rnaseq_pipeline_flowchart.pdf\n")
cat("- rnaseq_pipeline_flowchart.png\n") 