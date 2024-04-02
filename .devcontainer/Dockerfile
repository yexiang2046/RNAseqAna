FROM bioconductor/bioconductor_docker:RELEASE_3_17

# install R packages needed for VSCode interaction and package management
RUN R -e "install.packages(c('languageserver', 'renv', 'remotes', 'httpgd'))"

# install R packages for data science
RUN R -e 'install.packages(c("devtools", "viridis", "tidyverse", "hdf5r"))'


# install Bioinformatic tools
RUN R -e 'install.packages(c("pheatmap", "enrichR", "BiocManager", "ggfortify", "factoextra"))'
RUN R -e 'BiocManager::install(c("Rsubread", "edgeR", "DiffBind", "ChIPseeker", "ChIPpeakAnno", "ChIPQC", "limma", "msigdbr", "clusterProfiler", "fgsea", "org.Mm.eg.db", "org.Hs.eg.db", "Organism.dplyr", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Mmusculus.UCSC.mm39.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GSVA", "Rsubread", "profileplyr", "BSgenome.Mmusculus.UCSC.mm10"))'



# install radian via python and pip3
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends python3-setuptools \
    build-essential \
    bedtools \
    wget \
    bzip2 \
    ca-certificates \
    git \
    gcc \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libbz2-dev \
    libffi-dev \
    libgdbm-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss3-dev \
    libreadline-dev \
    libssl-dev \
    liblzma-dev \
    make \
    zlib1g-dev \
    rna-star \
    htop

RUN pip3 install radian 

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

RUN cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 &&   tar --bzip2 -xf samtools-1.9.tar.bz2

WORKDIR /tmp/samtools-1.9
RUN ./configure --enable-plugins --prefix=$SAMTOOLS_INSTALL_DIR &&   make all all-htslib &&   make install install-htslib
WORKDIR /
RUN ln -s $SAMTOOLS_INSTALL_DIR/bin/samtools /usr/bin/samtools &&   rm -rf /tmp/samtools-1.9
