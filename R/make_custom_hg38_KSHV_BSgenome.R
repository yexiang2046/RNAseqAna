# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")
library(BSgenome)

# Install other necessary packages
BiocManager::install("rtracklayer")

# Set the path to your FASTA file
fasta_file <- "path/to/your/fastafile.fasta"

# Metadata for GRCh38 and KSHV genome
organism <- "hg38_KSHV"
common_name <- "Human and Kaposi's sarcoma-associated herpesvirus"
provider <- "XiangYe"
provider_version <- "Version_of_the_data"
release_date <- "6/20/2024"
release_name <- "hg38_KSHV"
source_url <- "URL_to_the_source_data"

# Create a seed file for your genome
seed_file <- tempfile()
writeLines(con = seed_file, text = c(
  "Package: BSgenome.hg38.KSHV",
  "Title: Full genome sequences for hg38 and KSHV",
  "Description: This is a full genome sequence of the hg38 and KSHV genomes.",
  "Version: 1.0.0",
  "organism: human and Kaposi's sarcoma-associated herpesvirus",
  "common_name: CommonName",
  "provider: Provider",
  "provider_version: ProviderVersion",
  "release_date: ReleaseDate",
  "release_name: ReleaseName",
  "source_url: SourceURL",
  paste("organism_biocview: ", c("hg38", "KSHV"), sep="_"),
  "BSgenomeObjname: hg38_KSHV",
  "seqnames: c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'NC_009333.1')",  # Modify this line based on actual chromosome names in your FASTA
  "circ_seqs: character(0)",  # Update if there are circular chromosomes
  "single_sequences: character(0)",  # Update if there are unplaced scaffolds
  "masked_seqs: character(0)",
  "seqfile_name: GRCh38.primary_assembly.KSHV.genome.fa",
  "seqs_srcdir: ."
))

# Forge the BSgenome data package
forgeBSgenomeDataPkg(seed_file)

# This will create a directory in your current working directory that contains the new BSgenome package
# You may then install it using devtools or as you would a local R package
