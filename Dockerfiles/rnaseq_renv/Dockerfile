# Dockerfile for RNA-seq analysis with R/Bioconductor
# Includes edgeR, clusterProfiler, functional enrichment tools
# Version: v1.0.2 - Added ORA HALLMARK, geneset overlap filtering, gene deduplication

FROM rocker/r-ver:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN R -e "install.packages(c('BiocManager', 'optparse', 'ggplot2', 'pheatmap', 'data.table'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c(\
    'edgeR', \
    'clusterProfiler', \
    'enrichplot', \
    'fgsea', \
    'msigdbr', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'ReactomePA', \
    'factoextra'\
), update=TRUE, ask=FALSE)"

# Create directory for scripts
RUN mkdir -p /opt/bin

# Copy R scripts with new features
COPY bin/edger.r /opt/bin/edger.r
COPY bin/functional_analysis.r /opt/bin/functional_analysis.r

# Make scripts executable
RUN chmod +x /opt/bin/*.r

# Set working directory
WORKDIR /data

# Add scripts to PATH
ENV PATH="/opt/bin:${PATH}"

# Default command
CMD ["/bin/bash"]
