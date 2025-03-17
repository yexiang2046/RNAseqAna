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
             help = "Species for GO analysis (human or mouse)")
)

# Parse command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Function to extract gene information from GTF
extract_gene_info <- function(gtf_file) {
  # Read GTF file
  gtf_lines <- readLines(gtf_file)
  # Filter for gene entries
  gene_lines <- gtf_lines[grep('gene_id', gtf_lines) & grep('\tgene\t', gtf_lines)]
  
  # Extract gene_id and gene_name
  gene_ids <- gsub('.*gene_id "(.*?)".*', '\\1', gene_lines)
  gene_names <- gsub('.*gene_name "(.*?)".*', '\\1', gene_lines)
  
  # Create data frame
  gene_info <- data.frame(
    gene_id = gene_ids,
    gene_name = gene_names,
    stringsAsFactors = FALSE
  )
  
  return(gene_info)
}

# Load count data and metadata
counts <- read.delim(opt$counts, header = TRUE, skip = 1)
colnames(counts) <- sapply(colnames(counts), function(x){strsplit(x, "_")}[[1]][1])
colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- gsub("^X", "", colnames(counts))
# Extract gene expression counts (columns 7 onwards) 
counts_matrix <- counts[,7:ncol(counts)]

# Get gene annotations from GTF
gene_info <- extract_gene_info(opt$gtf)
# Match gene IDs from counts with GTF annotations
row.names(counts_matrix) <- counts$Geneid  # Assuming Geneid is the column with Ensembl IDs
gene_annotations <- data.frame(
  gene_id = counts$Geneid,
  gene_name = gene_info$gene_name[match(counts$Geneid, gene_info$gene_id)]
)

metaData <- read.table(opt$metadata, header = TRUE)

# order counts_matrix to metaData sample order with SampleId
counts_matrix <- counts_matrix[, metaData$SampleId]



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
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit the model
fit <- glmQLFit(y, design)

# Perform pairwise comparisons
group_levels <- levels(group)
comparison_results <- list()

for (i in 1:(length(group_levels) - 1)) {
  for (j in (i + 1):length(group_levels)) {
    contrast_name <- paste(group_levels[i], "vs", group_levels[j], sep = "_")
    contrast <- makeContrasts(contrasts = paste(group_levels[i], "-", group_levels[j]), levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    results <- topTags(qlf, n = Inf)
    # add annotation to results
    results$table$gene_id <- rownames(results$table)
    results$table$gene_name <- filtered_gene_annotations$gene_name[match(results$table$gene_id, filtered_gene_annotations$gene_id)]
    comparison_results[[contrast_name]] <- results
    write.csv(results, file = file.path(opt$output, paste0("DEG_", contrast_name, ".csv")), row.names = FALSE)
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

# Bar graph of differentially expressed genes for each comparison
for (contrast_name in names(comparison_results)) {
  de_genes <- comparison_results[[contrast_name]]$table
  de_genes$threshold <- as.factor(
    ifelse(de_genes$FDR < 0.05 & abs(de_genes$logFC) > 1,
           "Significant", 
           "Not Significant")
  )
  de_counts <- table(de_genes$threshold)
  bar_plot <- ggplot(as.data.frame(de_counts), 
                     aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Number of Differentially Expressed Genes:", 
                      contrast_name),
         x = "Group",
         y = "Count")
  ggsave(filename = file.path(opt$output, 
                             paste0("DEG_barplot_", contrast_name, ".png")),
         plot = bar_plot)
}

# Create a heatmap of the differentially expressed genes
# Get all the differentially expressed genes by FDR < 0.05
all_de_genes <- do.call(rbind, lapply(comparison_results, function(x) {
  x$table[x$table$FDR < 0.05 & abs(x$table$logFC) > 1, ]
}))

# remove duplicates
all_de_genes <- all_de_genes[!duplicated(all_de_genes$gene_id), ]$gene_id

# get the normalized counts for the differentially expressed genes
de_counts <- normalized_counts_with_annotations[all_de_genes, ]

# create a heatmap of the differentially expressed genes
heatmap_data <- de_counts[, -c(1:2)]
heatmap_data <- as.matrix(heatmap_data)
# scale the data by row
heatmap_data <- t(scale(t(heatmap_data)))

# annotation column
annotation_col <- data.frame(
  cell_type = metaData$Celltype,
  temperature = metaData$Temp,
  genotype = metaData$genotype
)
rownames(annotation_col) <- metaData$SampleId

# create a heatmap of the differentially expressed genes
heatmap_plot <- pheatmap(heatmap_data, 
                        #color = colorRampPalette(c("blue", "white", "red"))(100),
                        annotation_col = annotation_col,
                        show_rownames = FALSE,
                        cluster_rows = TRUE,
                        cluster_cols = TRUE)
ggsave(filename = file.path(opt$output, "DEG_heatmap.png"), plot = heatmap_plot)

# After the heatmap code, add k-means clustering analysis

# Function to perform k-means clustering and create plots
perform_kmeans_analysis <- function(data, n_clusters, label, output_dir) {
    # Scale the data
    scaled_data <- t(scale(t(data)))
    
    # Determine optimal number of clusters using elbow method
    wss <- sapply(1:15, function(k) {
        kmeans(scaled_data, centers=k, nstart=25)$tot.withinss
    })
    
    # Plot elbow curve
    elbow_plot <- ggplot(data.frame(k=1:15, wss=wss), aes(x=k, y=wss)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(title=paste("Elbow Method for", label),
             x="Number of Clusters (k)",
             y="Total Within Sum of Squares")
    ggsave(filename=file.path(output_dir, paste0("kmeans_elbow_", label, ".png")), 
           plot=elbow_plot)
    
    # Perform k-means clustering
    set.seed(42)
    km <- kmeans(scaled_data, centers=n_clusters, nstart=25)
    
    # Add cluster information to the data
    clustered_data <- data.frame(
        gene_id = rownames(data),
        cluster = km$cluster
    )
    
    # Create heatmap with cluster annotation
    annotation_row <- data.frame(
        Cluster = factor(km$cluster)
    )
    rownames(annotation_row) <- rownames(data)
    
    # Generate heatmap with clusters
    cluster_heatmap <- pheatmap(scaled_data,
                               annotation_col = annotation_col,
                               annotation_row = annotation_row,
                               show_rownames = FALSE,
                               cluster_rows = TRUE,
                               cluster_cols = TRUE,
                               main = paste("K-means Clustering (k=", n_clusters, ")", label))
    
    ggsave(filename=file.path(output_dir, paste0("kmeans_heatmap_", label, ".png")), 
           plot=cluster_heatmap)
    
    # Export cluster assignments
    write.csv(clustered_data, 
              file=file.path(output_dir, paste0("kmeans_clusters_", label, ".csv")),
              row.names=FALSE)
    
    # Create cluster profile plots
    cluster_means <- aggregate(scaled_data, 
                             by=list(Cluster=km$cluster), 
                             FUN=mean)
    
    # Melt the data for plotting
    cluster_means_long <- reshape2::melt(cluster_means, 
                                       id.vars="Cluster",
                                       variable.name="Sample",
                                       value.name="Expression")
    
    # Create profile plot
    profile_plot <- ggplot(cluster_means_long, 
                          aes(x=Sample, y=Expression, group=Cluster, color=factor(Cluster))) +
        geom_line() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title=paste("Cluster Profiles -", label),
             x="Samples",
             y="Scaled Expression",
             color="Cluster")
    
    ggsave(filename=file.path(output_dir, paste0("kmeans_profiles_", label, ".png")), 
           plot=profile_plot)
    
    return(km)
}

# Create directory for clustering results
clustering_dir <- file.path(opt$output, "clustering")
dir.create(clustering_dir, showWarnings = FALSE)

# Perform k-means clustering on all genes
all_genes_data <- normalized_counts_with_annotations[, -c(1:2)]
kmeans_all <- perform_kmeans_analysis(all_genes_data, 
                                    n_clusters=6, 
                                    label="all_genes",
                                    clustering_dir)

# Perform k-means clustering on DEGs only
deg_data <- normalized_counts_with_annotations[all_de_genes, -c(1:2)]
kmeans_deg <- perform_kmeans_analysis(deg_data, 
                                    n_clusters=4, 
                                    label="DEGs",
                                    clustering_dir)

# Modify the GO analysis section to include gene names and trim Ensembl IDs
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    # Load appropriate organism package based on species
    if (opt$species == "human") {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
    } else if (opt$species == "mouse") {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
    } else {
        stop("Unsupported species. Please use 'human' or 'mouse'")
    }
    
    # Create directory for GO analysis
    go_dir <- file.path(clustering_dir, "GO_analysis")
    dir.create(go_dir, showWarnings = FALSE)
    
    # Function to trim Ensembl ID version numbers
    trim_ensembl <- function(ids) {
        gsub("\\.[0-9]+$", "", ids)
    }
    
    # Function to perform GO analysis for clusters
    perform_go_analysis <- function(kmeans_result, gene_data, label) {
        for (i in 1:max(kmeans_result$cluster)) {
            # Get genes in current cluster
            cluster_genes <- rownames(gene_data)[kmeans_result$cluster == i]
            
            # Get gene symbols and descriptions for the cluster
            gene_symbols <- gene_annotations$gene_name[match(cluster_genes, gene_annotations$gene_id)]
            
            # Create and export annotated cluster gene list with expression values
            cluster_df <- data.frame(
                gene_id = cluster_genes,
                gene_name = gene_symbols,
                cluster = i,
                # Add expression values for each sample
                gene_data[cluster_genes, ]
            )
            
            # Export detailed cluster gene list
            write.csv(cluster_df,
                     file = file.path(go_dir, 
                                    paste0("cluster_", i, "_genes_", label, "_with_expression.csv")),
                     row.names = FALSE)
            
            # Trim Ensembl IDs for GO analysis
            trimmed_genes <- trim_ensembl(cluster_genes)
            
            # Perform GO analysis with trimmed IDs
            ego <- enrichGO(gene = trimmed_genes,
                          OrgDb = org_db,
                          keyType = "ENSEMBL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)
            
            # Save results
            if (nrow(ego) > 0) {
                # Add gene symbols to GO results
                ego_df <- as.data.frame(ego)
                genes_in_terms <- strsplit(ego_df$geneID, "/")
                ego_df$gene_symbols <- sapply(genes_in_terms, function(x) {
                    symbols <- gene_annotations$gene_name[match(x, trim_ensembl(gene_annotations$gene_id))]
                    paste(na.omit(symbols), collapse = "/")
                })
                
                write.csv(ego_df, 
                         file = file.path(go_dir, 
                                        paste0("GO_enrichment_", 
                                              label, "_cluster", i, ".csv")),
                         row.names = FALSE)
                
                # Create dot plot
                dot_plot <- dotplot(ego, showCategory = 20) +
                    ggtitle(paste("GO Enrichment -", label, "Cluster", i))
                ggsave(filename = file.path(go_dir, 
                                          paste0("GO_dotplot_", 
                                                label, "_cluster", i, ".png")),
                       plot = dot_plot)
            }
            
            # Create a simplified gene list with just IDs and symbols
            write.csv(data.frame(
                gene_id = cluster_genes,
                gene_name = gene_symbols,
                cluster = i
            ), file = file.path(go_dir, 
                              paste0("cluster_", i, "_genes_", label, "_simple.csv")),
            row.names = FALSE)
        }
    }
    
    # Perform GO analysis for both clustering results
    perform_go_analysis(kmeans_all, all_genes_data, "all_genes")
    perform_go_analysis(kmeans_deg, deg_data, "DEGs")
}





