# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")
library(BSgenome)

# Install other necessary packages
BiocManager::install("rtracklayer")

# Set the path to your FASTA file
fasta_file <- "GRCh38.primary_assembly.KSHV.genome.fa"

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
  #"provider_version: ProviderVersion",
  "release_date: ReleaseDate",
  #"release_name: ReleaseName",
  #"source_url: SourceURL",
  paste("organism_biocview: ", c("hg38", "KSHV"), sep="_"),
  "BSgenomeObjname: hg38_KSHV",
  "seqnames: c('chr10',  'chr12',  'chr15',  'chr18',  'chr21',  'chr3',  'chr6',  'chr9',   'chrY', 'chr11',  'chr13',  'chr16',  'chr19',  'chr22',  'chr4',  'chr7',  'chrM',  'NC_009333.1', 'chr1', 'chr14', 'chr17', 'chr20', 'chr2', 'chr5', 'chr8', 'chrX')",  # Modify this line based on actual chromosome names in your FASTA
  #"circ_seqs: character(0)",  # Update if there are circular chromosomes
  #"single_sequences: character(0)",  # Update if there are unplaced scaffolds
  #"masked_seqs: character(0)",
  "genome: GRCh38.p13",
  "seqfile_name: c('chr10.fasta',  'chr12.fasta',  'chr15.fasta',  'chr18.fasta',  'chr21.fasta',  'chr3.fasta',  'chr6.fasta',  'chr9.fasta',   'chrY.fasta', 'chr11.fasta',  'chr13.fasta',  'chr16.fasta',  'chr19.fasta',  'chr22.fasta',  'chr4.fasta',  'chr7.fasta',  'chrM.fasta',  'NC_009333.1.fasta', 'chr1.fasta', 'chr14.fasta', 'chr17.fasta', 'chr20.fasta', 'chr2.fasta', 'chr5.fasta', 'chr8.fasta', 'chrX.fasta')",
  "seqs_srcdir: chr_fasta"
))

# Forge the BSgenome data package
forgeBSgenomeDataPkg(seed_file)

# This will create a directory in your current working directory that contains the new BSgenome package
# You may then install it using devtools or as you would a local R package
