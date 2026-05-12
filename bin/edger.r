#!/usr/bin/env Rscript

# Load necessary libraries
library(edgeR)
library(optparse)
library(ggplot2)
library(factoextra)
library(patchwork)

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to counts file"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to metadata file"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for results"),
  make_option(c("-g", "--gtf"), type = "character", help = "Path to GTF annotation file"),
  make_option(c("-s", "--species"), type = "character", default = "human", help = "Species for GO analysis (human or mouse)"),
  make_option(c("-p", "--paired"), type = "character", default = NULL, help = "Column name in metadata for paired/blocking factor (e.g. 'patient_id')"),
  make_option(c("-f", "--factors"), type = "character", default = NULL, help = "Comma-separated list of metadata columns to combine with 'group' (e.g. 'time')")
)

# Parse command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Function to extract gene information from GTF
extract_gene_info <- function(gtf_file) {
  gtf_lines <- readLines(gtf_file)
  gene_lines <- gtf_lines[grepl('gene_id', gtf_lines) & grepl('\tgene\t', gtf_lines)]
  
  gene_ids   <- gsub('.*gene_id "(.*?)".*', '\\1', gene_lines)
  gene_names <- gsub('.*gene_name "(.*?)".*', '\\1', gene_lines)
  gene_types <- gsub('.*gene_type "(.*?)".*', '\\1', gene_lines)
  
  data.frame(gene_id = gene_ids, gene_name = gene_names, gene_type = gene_types, stringsAsFactors = FALSE)
}

# === LOAD COUNTS ===
counts <- read.delim(opt$counts, header = TRUE, skip = 1, check.names = FALSE)

cat("=== RAW column names from counts.txt (first 6) ===\n")
print(head(colnames(counts), 6))

col_names <- basename(colnames(counts))
col_names <- gsub("Aligned\\.sortedByCoord\\.out\\.bam$", "", col_names)
col_names <- gsub("\\.sortedByCoord\\.out\\.bam$", "", col_names)
col_names <- gsub("\\.bam$", "", col_names)
col_names <- gsub("^X", "", col_names)

colnames(counts) <- col_names
counts_matrix <- counts[, 7:ncol(counts)]

cat("\n=== CLEANED sample names (first 6) ===\n")
print(head(colnames(counts_matrix), 6))

# Gene annotations
gene_info <- extract_gene_info(opt$gtf)
row.names(counts_matrix) <- counts$Geneid

gene_annotations <- data.frame(
  gene_id   = counts$Geneid,
  gene_name = gene_info$gene_name[match(counts$Geneid, gene_info$gene_id)],
  gene_type = gene_info$gene_type[match(counts$Geneid, gene_info$gene_id)]
)

# Remove duplicate gene names
n_before <- nrow(gene_annotations)
dup_gene_names <- unique(gene_annotations$gene_name[duplicated(gene_annotations$gene_name) | duplicated(gene_annotations$gene_name, fromLast = TRUE)])
dup_gene_names <- dup_gene_names[!is.na(dup_gene_names)]

if (length(dup_gene_names) > 0) {
  cat("Found", length(dup_gene_names), "duplicated gene names\n")
  total_expr <- rowSums(counts_matrix)
  genes_to_keep <- rep(TRUE, nrow(gene_annotations))
  for (dup_name in dup_gene_names) {
    dup_idx <- which(gene_annotations$gene_name == dup_name)
    if (length(dup_idx) > 1) {
      max_idx <- dup_idx[which.max(total_expr[dup_idx])]
      genes_to_keep[dup_idx] <- FALSE
      genes_to_keep[max_idx] <- TRUE
    }
  }
  counts_matrix    <- counts_matrix[genes_to_keep, ]
  gene_annotations <- gene_annotations[genes_to_keep, ]
  cat("Removed", n_before - nrow(gene_annotations), "duplicate gene entries\n")
}

# Load metadata
metaData <- read.table(opt$metadata, header = TRUE)
metaData$group <- factor(metaData$group)

# === BUILD COMBINED GROUP AS NEW MAIN GROUP ===
additional_factors <- character(0)
if (!is.null(opt$factors) && opt$factors != "") {
  additional_factors <- trimws(strsplit(opt$factors, ",")[[1]])
  missing <- setdiff(additional_factors, colnames(metaData))
  if (length(missing) > 0) {
    stop(paste("❌ ERROR: The following factors were not found in metadata.txt:", paste(missing, collapse = ", ")))
  }
  for (col in additional_factors) {
    metaData[[col]] <- factor(metaData[[col]])
  }
  cat("✅ Factors combined with group:", paste(additional_factors, collapse = ", "), "\n")
}

combined_vars <- c("group", additional_factors)
metaData$combined_group <- interaction(metaData[, combined_vars, drop = FALSE], sep = "_", drop = TRUE, lex.order = TRUE)
metaData$combined_group <- factor(metaData$combined_group)

cat("✅ Using COMBINED_GROUP as the new main group (", nlevels(metaData$combined_group), "levels):\n")
print(levels(metaData$combined_group))

# Order counts
missing_samples <- setdiff(metaData$SampleId, colnames(counts_matrix))
if (length(missing_samples) > 0) {
  cat("\n❌ ERROR: These SampleIds from metadata.txt are MISSING in the counts file:\n")
  print(missing_samples)
  stop("Sample name mismatch.")
}

counts_matrix <- counts_matrix[, metaData$SampleId]

# Create DGEList using combined_group
y <- DGEList(counts = counts_matrix, group = metaData$combined_group)

# Filter & Normalize
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Export normalized counts
normalized_counts <- cpm(y)
filtered_gene_annotations <- gene_annotations[keep, ]
row.names(normalized_counts) <- row.names(y)
row.names(filtered_gene_annotations) <- filtered_gene_annotations$gene_id
normalized_counts_with_annotations <- cbind(
  filtered_gene_annotations[row.names(normalized_counts), ],
  normalized_counts
)
write.csv(normalized_counts_with_annotations, file.path(opt$output, "normalized_counts_CPM.csv"), row.names = FALSE)

# === DESIGN MATRIX ===
if (!is.null(opt$paired)) {
  if (!(opt$paired %in% colnames(metaData))) stop(paste("❌ Paired column", opt$paired, "not found"))
  metaData[[opt$paired]] <- factor(metaData[[opt$paired]])
  design <- model.matrix(as.formula(paste0("~0 + ", opt$paired, " + combined_group")), data = metaData)
  cat("✅ PAIRED design using combined_group\n")
} else {
  design <- model.matrix(~0 + combined_group, data = metaData)
  cat("✅ UNPAIRED design using combined_group\n")
}

cat("Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")

y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# === ALL PAIRWISE COMPARISONS ===
cat("\n=== Performing pairwise comparisons on combined_group ===\n")
comb_levels <- levels(metaData$combined_group)
for (i in 1:(length(comb_levels)-1)) {
  for (j in (i+1):length(comb_levels)) {
    contrast_name <- paste(comb_levels[j], "vs", comb_levels[i], sep = "_")
    
    col_j <- paste0("combined_group", comb_levels[j])
    col_i <- paste0("combined_group", comb_levels[i])
    
    if (col_j %in% colnames(design) && col_i %in% colnames(design)) {
      contrast_str <- paste0("`", col_j, "` - `", col_i, "`")
    } else if (col_j %in% colnames(design)) {
      contrast_str <- paste0("`", col_j, "`")
    } else if (col_i %in% colnames(design)) {
      contrast_str <- paste0("-`", col_i, "`")
    } else {
      next
    }
    
    cat("  ", contrast_name, "→", contrast_str, "\n")
    contrast <- makeContrasts(contrasts = contrast_str, levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    results <- topTags(qlf, n = Inf)
    
    results$table$gene_id   <- rownames(results$table)
    results$table$gene_name <- filtered_gene_annotations$gene_name[match(results$table$gene_id, filtered_gene_annotations$gene_id)]
    results$table$gene_type <- filtered_gene_annotations$gene_type[match(results$table$gene_id, filtered_gene_annotations$gene_id)]
    
    write.csv(results$table, file.path(opt$output, paste0("DEG_", contrast_name, ".csv")), row.names = FALSE)
  }
}

# === PCA PLOTS - MULTIPLE PANELS WITH COMBINED GROUP ===
pca_data <- prcomp(t(cpm(y, log = TRUE)), scale. = TRUE)
var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100

plot_pca_pair <- function(ax1, ax2, color_by, title) {
  fviz_pca_ind(pca_data, axes = c(ax1, ax2), label = "none",
               habillage = color_by, addEllipses = TRUE, ellipse.level = 0.95,
               mean.point = FALSE) +
    labs(x = sprintf("PC%d (%.1f%%)", ax1, var_explained[ax1]),
         y = sprintf("PC%d (%.1f%%)", ax2, var_explained[ax2])) +
    ggtitle(title) + theme_classic()
}

# 1. Multi-panel PCA colored by ORIGINAL group (kept for reference)
p12 <- plot_pca_pair(1, 2, metaData$group, "PC1 vs PC2 - by original group")
p13 <- plot_pca_pair(1, 3, metaData$group, "PC1 vs PC3 - by original group")
p14 <- plot_pca_pair(1, 4, metaData$group, "PC1 vs PC4 - by original group")
p23 <- plot_pca_pair(2, 3, metaData$group, "PC2 vs PC3 - by original group")
p24 <- plot_pca_pair(2, 4, metaData$group, "PC2 vs PC4 - by original group")
p34 <- plot_pca_pair(3, 4, metaData$group, "PC3 vs PC4 - by original group")

pca_multi_orig <- (p12 + p13 + p14) / (p23 + p24 + p34) +
  plot_annotation(title = "PCA - Colored by Original Group", theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(file.path(opt$output, "PCA_multi_plot_original_group.png"), pca_multi_orig, width = 15, height = 10, dpi = 300)

# 2. NEW: Multi-panel PCA colored by COMBINED GROUP (what you asked for)
p12_comb <- plot_pca_pair(1, 2, metaData$combined_group, "PC1 vs PC2 - by combined_group")
p13_comb <- plot_pca_pair(1, 3, metaData$combined_group, "PC1 vs PC3 - by combined_group")
p14_comb <- plot_pca_pair(1, 4, metaData$combined_group, "PC1 vs PC4 - by combined_group")
p23_comb <- plot_pca_pair(2, 3, metaData$combined_group, "PC2 vs PC3 - by combined_group")
p24_comb <- plot_pca_pair(2, 4, metaData$combined_group, "PC2 vs PC4 - by combined_group")
p34_comb <- plot_pca_pair(3, 4, metaData$combined_group, "PC3 vs PC4 - by combined_group")

pca_multi_comb <- (p12_comb + p13_comb + p14_comb) / (p23_comb + p24_comb + p34_comb) +
  plot_annotation(title = "PCA - Colored by Combined Group (new main group)", theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(file.path(opt$output, "PCA_multi_plot_combined_group.png"), pca_multi_comb, width = 15, height = 10, dpi = 300)

# 3. Single high-res PCA colored by combined_group
p_comb <- fviz_pca_ind(pca_data, axes = c(1, 2), label = "none",
                       habillage = metaData$combined_group, addEllipses = TRUE,
                       mean.point = FALSE) +
          ggtitle("PCA colored by combined_group (new main group)") + theme_classic()
ggsave(file.path(opt$output, "PCA_by_combined_group.png"), p_comb, width = 9, height = 7, dpi = 300)