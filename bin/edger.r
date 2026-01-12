#!/usr/bin/env Rscript

# Load necessary libraries
library(edgeR)
library(optparse)
library(ggplot2)
library(factoextra)
library(pheatmap)
library(clusterProfiler)

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to counts file"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to metadata file"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for results"),
  make_option(c("-g", "--gtf"), type = "character", help = "Path to GTF annotation file"),
  make_option(c("-s", "--species"), type = "character", default = "human",
             help = "Species for GO analysis (human or mouse)"),
  make_option(c("-p", "--paired"), type = "character", default = NULL,
             help = "Column name in metadata for paired samples (e.g., 'patient_id', 'subject')")
)

# Parse command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Function to extract gene information from GTF
extract_gene_info <- function(gtf_file) {
  # Read GTF file
  gtf_lines <- readLines(gtf_file)
  # Filter for gene entries
  gene_lines <- gtf_lines[grep('gene_id', gtf_lines) & grep('\tgene\t', gtf_lines)]
  
  # Extract gene_id, gene_name, and gene_type
  gene_ids <- gsub('.*gene_id "(.*?)".*', '\\1', gene_lines)
  gene_names <- gsub('.*gene_name "(.*?)".*', '\\1', gene_lines)
  gene_types <- gsub('.*gene_type "(.*?)".*', '\\1', gene_lines)
  
  # Create data frame
  gene_info <- data.frame(
    gene_id = gene_ids,
    gene_name = gene_names,
    gene_type = gene_types,
    stringsAsFactors = FALSE
  )
  
  return(gene_info)
}

# Load count data and metadata
counts <- read.delim(opt$counts, header = TRUE, skip = 1)
# R converts column names: X14197.XZ.0003_S1_L005Aligned... -> X14197.XZ.0003_S1_L005Aligned...
# Convert X14197.XZ.0003_S1_L005Aligned.sortedByCoord.out.bam -> 14197-XZ-3
colnames(counts) <- gsub("^X(\\d+)\\.(\\w+)\\.0*(\\d+)_.*", "\\1-\\2-\\3", colnames(counts))
# Extract gene expression counts (columns 7 onwards) 
counts_matrix <- counts[,7:ncol(counts)]

# Get gene annotations from GTF
gene_info <- extract_gene_info(opt$gtf)
# Match gene IDs from counts with GTF annotations
row.names(counts_matrix) <- counts$Geneid  # Assuming Geneid is the column with Ensembl IDs
gene_annotations <- data.frame(
  gene_id = counts$Geneid,
  gene_name = gene_info$gene_name[match(counts$Geneid, gene_info$gene_id)],
  gene_type = gene_info$gene_type[match(counts$Geneid, gene_info$gene_id)]
)

metaData <- read.table(opt$metadata, header = TRUE)

# order counts_matrix to metaData sample order with SampleId
counts_matrix <- counts_matrix[, metaData$SampleId]

# Read metadata and remove SampleId column
metaData_for_annotation <- metaData[, !colnames(metaData) %in% "SampleId"]


# Create annotation column using all remaining columns for heatmap
annotation_col <- as.data.frame(metaData_for_annotation)
rownames(annotation_col) <- metaData$SampleId


# Create DGEList object
group <- factor(metaData$group)
y <- DGEList(counts = counts_matrix, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Export normalized counts (CPM)
normalized_counts <- cpm(y)
# Add gene annotations back to the normalized counts
filtered_gene_annotations <- gene_annotations[keep, ]
# Ensure gene IDs match by setting row names
row.names(normalized_counts) <- row.names(y)
row.names(filtered_gene_annotations) <- filtered_gene_annotations$gene_id
# Verify and combine in correct order
normalized_counts_with_annotations <- cbind(
  filtered_gene_annotations[row.names(normalized_counts), ],
  normalized_counts
)
write.csv(normalized_counts_with_annotations,
          file = file.path(opt$output, "normalized_counts_CPM.csv"),
          row.names = FALSE)

# Create design matrix
# Check if paired analysis is requested
if (!is.null(opt$paired)) {
  # Verify the paired column exists in metadata
  if (!(opt$paired %in% colnames(metaData))) {
    stop(paste("Paired column", opt$paired, "not found in metadata"))
  }

  # Create paired factor
  pair <- factor(metaData[[opt$paired]])

  # Create design matrix with paired samples (use ~0 to remove intercept for consistent column naming)
  design <- model.matrix(~0 + pair + group)

  cat("Performing paired analysis using column:", opt$paired, "\n")
  cat("Design matrix formula: ~0 + pair + group\n")
  cat("Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")
} else {
  # Unpaired design (original behavior)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  cat("Performing unpaired analysis\n")
}

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit the model
fit <- glmQLFit(y, design)

# Perform pairwise comparisons
group_levels <- levels(group)
comparison_results <- list()

if (!is.null(opt$paired)) {
  # Paired design: perform all pairwise comparisons using contrasts
  # The design matrix includes pair effects, so contrasts control for pairing
  # Note: reference group level won't have a column, so we need to handle that
  
  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {
      contrast_name <- paste(group_levels[j], "vs", group_levels[i], sep = "_")
      
      # Build contrast string: if group is in design matrix use it, otherwise use negative of others
      group_j_col <- paste0("group", group_levels[j])
      group_i_col <- paste0("group", group_levels[i])
      
      cat("Comparing:", contrast_name, "\n")
      cat("  group_j_col:", group_j_col, "exists:", (group_j_col %in% colnames(design)), "\n")
      cat("  group_i_col:", group_i_col, "exists:", (group_i_col %in% colnames(design)), "\n")
      
      # Check which columns exist in design matrix
      if (group_j_col %in% colnames(design) & group_i_col %in% colnames(design)) {
        contrast_str <- paste0("`", group_j_col, "` - `", group_i_col, "`")
      } else if (group_j_col %in% colnames(design)) {
        # group_i is reference, so contrast is just the group_j coefficient
        contrast_str <- paste0("`", group_j_col, "`")
      } else if (group_i_col %in% colnames(design)) {
        # group_j is reference, so contrast is -group_i
        contrast_str <- paste0("-`", group_i_col, "`")
      } else {
        # Both are reference level, skip
        cat("  Skipping: both are reference level\n")
        next
      }
      
      cat("  contrast_str:", contrast_str, "\n")
      contrast <- makeContrasts(contrasts = contrast_str, levels = design)
      qlf <- glmQLFTest(fit, contrast = contrast)
      results <- topTags(qlf, n = Inf)

      # Add annotation to results
      results$table$gene_id <- rownames(results$table)
      results$table$gene_name <- filtered_gene_annotations$gene_name[match(results$table$gene_id, filtered_gene_annotations$gene_id)]
      results$table$gene_type <- filtered_gene_annotations$gene_type[match(results$table$gene_id, filtered_gene_annotations$gene_id)]

      comparison_results[[contrast_name]] <- results
      write.csv(results, file = file.path(opt$output, paste0("DEG_", contrast_name, ".csv")), row.names = FALSE)
    }
  }
} else {
  # Unpaired design: use contrasts for all pairwise comparisons
  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {
      contrast_name <- paste(group_levels[j], "vs", group_levels[i], sep = "_")
      contrast_str <- paste0("`", group_levels[j], "` - `", group_levels[i], "`")
      contrast <- makeContrasts(contrasts = contrast_str, levels = design)
      qlf <- glmQLFTest(fit, contrast = contrast)
      results <- topTags(qlf, n = Inf)

      # Add annotation to results
      results$table$gene_id <- rownames(results$table)
      results$table$gene_name <- filtered_gene_annotations$gene_name[match(results$table$gene_id, filtered_gene_annotations$gene_id)]
      results$table$gene_type <- filtered_gene_annotations$gene_type[match(results$table$gene_id, filtered_gene_annotations$gene_id)]

      comparison_results[[contrast_name]] <- results
      write.csv(results, file = file.path(opt$output, paste0("DEG_", contrast_name, ".csv")), row.names = FALSE)
    }
  }
}

# Plot PCA
pca_data <- prcomp(t(cpm(y, log = TRUE)), scale. = TRUE)
pca_plot <- fviz_pca_ind(pca_data, 
                         label = "none", 
                         habillage = group, 
                         addEllipses = FALSE) +
  ggtitle("PCA of Samples") +
  theme_classic()
ggsave(filename = file.path(opt$output, "PCA_plot.png"), plot = pca_plot)
