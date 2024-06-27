# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")
library(BSgenome)


forgeBSgenomeDataPkg("BSgenome_seedfile.txt")


