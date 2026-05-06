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

# === LOAD COUNTS WITH ROBUST NAME CLEANING + DIAGNOSTICS ===
counts <- read.delim(opt$counts, header = TRUE, skip = 1, check.names = FALSE)

# Show raw names (so we can debug)
cat("=== RAW column names from counts.txt (first 6) ===\n")
print(head(colnames(counts), 6))

# Clean the column names properly
col_names <- basename(colnames(counts))                    # remove any folder path
col_names <- gsub("Aligned\\.sortedByCoord\\.out\\.bam$", "", col_names)
col_names <- gsub("\\.sortedByCoord\\.out\\.bam$", "", col_names)
col_names <- gsub("\\.bam$", "", col_names)
col_names <- gsub("^X", "", col_names)                    # remove R's occasional X prefix

colnames(counts) <- col_names
counts_matrix <- counts[, 7:ncol(counts)]

cat("\n=== CLEANED sample names (first 6) ===\n")
print(head(colnames(counts_matrix), 6))


# Get gene annotations from GTF
gene_info <- extract_gene_info(opt$gtf)
# Match gene IDs from counts with GTF annotations
row.names(counts_matrix) <- counts$Geneid  # Assuming Geneid is the column with Ensembl IDs
gene_annotations <- data.frame(
  gene_id = counts$Geneid,
  gene_name = gene_info$gene_name[match(counts$Geneid, gene_info$gene_id)],
  gene_type = gene_info$gene_type[match(counts$Geneid, gene_info$gene_id)]
)

# Remove duplicate gene names by keeping the one with highest total expression
# This prevents issues in downstream enrichment analyses
n_before <- nrow(gene_annotations)
dup_gene_names <- gene_annotations$gene_name[duplicated(gene_annotations$gene_name) |
                                              duplicated(gene_annotations$gene_name, fromLast = TRUE)]
dup_gene_names <- unique(dup_gene_names[!is.na(dup_gene_names)])

if (length(dup_gene_names) > 0) {
  cat("Found", length(dup_gene_names), "duplicated gene names\n")

  # Calculate total expression for each gene
  total_expr <- rowSums(counts_matrix)

  # For each duplicated gene name, keep only the one with highest expression
  genes_to_keep <- rep(TRUE, nrow(gene_annotations))
  for (dup_name in dup_gene_names) {
    dup_indices <- which(gene_annotations$gene_name == dup_name)
    if (length(dup_indices) > 1) {
      # Keep the one with highest total expression
      max_expr_idx <- dup_indices[which.max(total_expr[dup_indices])]
      genes_to_keep[dup_indices] <- FALSE
      genes_to_keep[max_expr_idx] <- TRUE
    }
  }

  # Filter counts_matrix and gene_annotations
  counts_matrix <- counts_matrix[genes_to_keep, ]
  gene_annotations <- gene_annotations[genes_to_keep, ]

  n_after <- nrow(gene_annotations)
  cat("Removed", n_before - n_after, "duplicate gene entries (kept highest expressed)\n")
}

metaData <- read.table(opt$metadata, header = TRUE)

# Safety check: do all metadata samples exist in the counts file?
missing <- setdiff(metaData$SampleId, colnames(counts_matrix))
if (length(missing) > 0) {
  cat("\n❌ ERROR: These SampleIds from metadata.txt are MISSING in the counts file:\n")
  print(missing)
  stop("Sample name mismatch between metadata and featureCounts output.")
} else {
  cat("\n✅ All metadata samples found in counts file.\n")
}

# Show metadata SampleIds for comparison
cat("\n=== Metadata SampleId (first 6) ===\n")
print(head(metaData$SampleId, 6))

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

# Plot PCA - All 6 pairwise combinations of first 4 PCs
pca_data <- prcomp(t(cpm(y, log = TRUE)), scale. = TRUE)

# Calculate percentage of variance explained by each PC
var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100

# Helper function to create a PCA plot for any pair of dimensions
plot_pca_pair <- function(ax1, ax2) {
  fviz_pca_ind(pca_data, 
               axes = c(ax1, ax2),
               label = "none", 
               habillage = group, 
               addEllipses = FALSE,
               mean.point = FALSE) +
    labs(x = sprintf("PC%d (%.1f%%)", ax1, var_explained[ax1]),
         y = sprintf("PC%d (%.1f%%)", ax2, var_explained[ax2])) +
    ggtitle(sprintf("PC%d vs PC%d", ax1, ax2)) +
    theme_classic()
}

# Generate all 6 pairwise plots (exactly C(4,2) = 6 combinations)
p12 <- plot_pca_pair(1, 2)
p13 <- plot_pca_pair(1, 3)
p14 <- plot_pca_pair(1, 4)
p23 <- plot_pca_pair(2, 3)
p24 <- plot_pca_pair(2, 4)
p34 <- plot_pca_pair(3, 4)

# Use patchwork for clean multi-panel layout (2 rows × 3 columns)
library(patchwork)

pca_multi <- (p12 + p13 + p14) / 
             (p23 + p24 + p34) +
  plot_annotation(title = "PCA of Samples - All Pairwise Plots (PC1 to PC4)",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

# Save the combined plot (larger size to fit 6 panels nicely)
ggsave(filename = file.path(opt$output, "PCA_multi_plot.png"), 
       plot = pca_multi, 
       width = 15, height = 10, dpi = 300)