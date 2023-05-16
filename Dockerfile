From ubuntu:18.04

LABEL For CAGEr analysis v1.0

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y upgrade \
    && apt-get install -y \
    wget
RUN apt-get update && apt-get -y upgrade \
    && apt-get install -y \
    build-essential \
    bzip2 \
    ca-certificates \
    dirmngr \
    gcc \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libglib2.0-0 \
    liblzma-dev \
    libncurses5-dev  \
    libncursesw5-dev \
    libsm6 \
    libxext6 \
    libxml2-dev \
    libxrender1 \
    make \
    software-properties-common \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* 

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/"
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev 

RUN apt-get update
   
RUN cd /usr/bin \
    && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && make

RUN cd .. \
    && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -vxjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make \
    && cd ..

RUN wget https://bioconductor.org/packages/3.16/bioc/src/contrib/BiocGenerics_0.44.0.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/BiocParallel_1.32.1.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/BSgenome_1.66.1.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/data.table_1.14.4.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/DelayedArray_0.24.0.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/DelayedMatrixStats_1.20.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/formula.tools_1.7.1.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/GenomeInfoDb_1.34.3.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/GenomicAlignments_1.34.0.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/GenomicRanges_1.50.1.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/ggplot2_3.4.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/gtools_3.9.3.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/IRanges_2.32.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/KernSmooth_2.23-20.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/memoise_2.0.1.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/plyr_1.8.8.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/Rsamtools_2.14.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/reshape2_1.4.4.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/rtracklayer_1.58.0.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/S4Vectors_0.36.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/som_0.3-5.1.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/stringdist_0.9.10.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/stringi_1.7.8.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/SummarizedExperiment_1.28.0.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/vegan_2.6-4.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/VGAM_1.1-7.tar.gz \
    && wget https://cran.r-project.org/src/contrib/futile.logger_1.4.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/snow_0.4-4.tar.gz \
    && wget https://cran.r-project.org/src/contrib/BH_1.78.0-0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/cpp11_0.4.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/lambda.r_1.2.4.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/futile.options_1.0.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/formatR_1.12.tar.gz \
    && wget https://cran.rstudio.com/src/contrib/matrixStats_0.62.0.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/MatrixGenerics_1.10.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/Rcpp_1.0.9.tar.gz \
    && wget https://www.bioconductor.org/packages/release/bioc/src/contrib/sparseMatrixStats_1.10.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/operator.tools_1.6.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/bitops_1.0-7.tar.gz \
    && wget https://cran.r-project.org/src/contrib/RCurl_1.98-1.9.tar.gz \
    && wget https://bioconductor.org/packages/release/data/annotation/src/contrib/GenomeInfoDbData_1.2.9.tar.gz \
    && wget https://cran.r-project.org/src/contrib/cli_3.4.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/glue_1.6.2.tar.gz \
    && wget https://cran.r-project.org/src/contrib/gtable_0.3.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/isoband_0.2.6.tar.gz \
    && wget https://cran.r-project.org/src/contrib/lifecycle_1.0.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/rlang_1.0.6.tar.gz \
    && wget https://cran.r-project.org/src/contrib/scales_1.2.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/tibble_3.1.8.tar.gz \
    && wget https://cran.r-project.org/src/contrib/vctrs_0.5.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/farver_2.1.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/labeling_0.4.2.tar.gz \
    && wget https://cran.r-project.org/src/contrib/munsell_0.5.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/R6_2.5.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/RColorBrewer_1.1-3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/viridisLite_0.4.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/colorspace_2.0-3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/fansi_1.0.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/magrittr_2.0.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/pillar_1.8.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/pkgconfig_2.0.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/utf8_1.2.2.tar.gz \
    && wget https://cran.r-project.org/src/contrib/withr_2.5.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/fastmap_1.1.0.tar.gz \ 
    && wget https://cran.r-project.org/src/contrib/cachem_1.0.6.tar.gz \
    && wget https://cran.r-project.org/src/contrib/stringr_1.4.1.tar.gz \
    && wget https://cran.r-project.org/src/contrib/permute_0.9-7.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/zlibbioc_1.44.0.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/XVector_0.38.0.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/Rhtslib_2.0.0.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.58.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/crayon_1.5.2.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/Biostrings_2.66.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/XML_3.99-0.12.tar.gz \
    && wget https://www.bioconductor.org/packages/release/bioc/src/contrib/BiocIO_1.8.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/rjson_0.2.21.tar.gz \
    && wget https://cran.r-project.org/src/contrib/yaml_2.3.6.tar.gz \
    && wget https://cran.r-project.org/src/contrib/restfulr_0.0.15.tar.gz \
    && wget https://bioconductor.org/packages/devel/bioc/src/contrib/BiocBaseUtils_1.1.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/tidyselect_1.2.0.tar.gz \
    && wget https://cran.r-project.org/src/contrib/generics_0.1.3.tar.gz \
    && wget https://cran.r-project.org/src/contrib/dplyr_1.0.10.tar.gz \
    && wget https://cran.r-project.org/src/contrib/ellipsis_0.3.2.tar.gz \
    && wget https://cran.r-project.org/src/contrib/purrr_0.3.5.tar.gz \
    && wget https://cran.r-project.org/src/contrib/tidyr_1.2.1.tar.gz \
    && wget https://bioconductor.org/packages/release/bioc/src/contrib/MultiAssayExperiment_1.24.0.tar.gz \
    && wget https://bioconductor.org/packages/3.16/bioc/src/contrib/CAGEr_2.4.0.tar.gz \
    && R CMD INSTALL BiocGenerics_0.44.0.tar.gz \
    formatR_1.12.tar.gz \
    lambda.r_1.2.4.tar.gz \
    futile.options_1.0.1.tar.gz \
    futile.logger_1.4.3.tar.gz \
    snow_0.4-4.tar.gz \
    BH_1.78.0-0.tar.gz \
    cpp11_0.4.3.tar.gz \
    BiocParallel_1.32.1.tar.gz \
    data.table_1.14.4.tar.gz \
    S4Vectors_0.36.0.tar.gz \
    IRanges_2.32.0.tar.gz \
    matrixStats_0.62.0.tar.gz \
    MatrixGenerics_1.10.0.tar.gz \
    DelayedArray_0.24.0.tar.gz \
    Rcpp_1.0.9.tar.gz \
    sparseMatrixStats_1.10.0.tar.gz \
    DelayedMatrixStats_1.20.0.tar.gz \
    operator.tools_1.6.3.tar.gz \
    formula.tools_1.7.1.tar.gz \
    bitops_1.0-7.tar.gz \
    RCurl_1.98-1.9.tar.gz \
    GenomeInfoDbData_1.2.9.tar.gz \
    GenomeInfoDb_1.34.3.tar.gz \
    ALL cli_3.4.1.tar.gz \
    glue_1.6.2.tar.gz \
    gtable_0.3.1.tar.gz \
    isoband_0.2.6.tar.gz \
    rlang_1.0.6.tar.gz \
    lifecycle_1.0.3.tar.gz \
    ALL farver_2.1.1.tar.gz \
    labeling_0.4.2.tar.gz \
    colorspace_2.0-3.tar.gz \
    munsell_0.5.0.tar.gz \
    R6_2.5.1.tar.gz \
    RColorBrewer_1.1-3.tar.gz \
    viridisLite_0.4.1.tar.gz \
    scales_1.2.1.tar.gz \
    vctrs_0.5.0.tar.gz \
    ALL fansi_1.0.3.tar.gz \
    magrittr_2.0.3.tar.gz \
    utf8_1.2.2.tar.gz \
    pillar_1.8.1.tar.gz \
    pkgconfig_2.0.3.tar.gz \
    tibble_3.1.8.tar.gz \
    withr_2.5.0.tar.gz \
    ggplot2_3.4.0.tar.gz \
    gtools_3.9.3.tar.gz \
    KernSmooth_2.23-20.tar.gz \
    fastmap_1.1.0.tar.gz \
    cachem_1.0.6.tar.gz \
    memoise_2.0.1.tar.gz \
    plyr_1.8.8.tar.gz \
    stringi_1.7.8.tar.gz \
    stringr_1.4.1.tar.gz \
    reshape2_1.4.4.tar.gz \
    som_0.3-5.1.tar.gz \
    stringdist_0.9.10.tar.gz \
    permute_0.9-7.tar.gz \
    vegan_2.6-4.tar.gz \
    VGAM_1.1-7.tar.gz \
    zlibbioc_1.44.0.tar.gz \
    XVector_0.38.0.tar.gz \
    Rhtslib_2.0.0.tar.gz \
    GenomicRanges_1.50.1.tar.gz \
    Biobase_2.58.0.tar.gz \
    SummarizedExperiment_1.28.0.tar.gz \
    crayon_1.5.2.tar.gz \
    Biostrings_2.66.0.tar.gz \
    Rsamtools_2.14.0.tar.gz \
    GenomicAlignments_1.34.0.tar.gz \
    XML_3.99-0.12.tar.gz \
    BiocIO_1.8.0.tar.gz \
    rjson_0.2.21.tar.gz \
    yaml_2.3.6.tar.gz \
    restfulr_0.0.15.tar.gz \
    rtracklayer_1.58.0.tar.gz \
    BSgenome_1.66.1.tar.gz \
    BiocBaseUtils_1.1.0.tar.gz \
    tidyselect_1.2.0.tar.gz \
    generics_0.1.3.tar.gz \
    dplyr_1.0.10.tar.gz \
    ellipsis_0.3.2.tar.gz \
    purrr_0.3.5.tar.gz \
    tidyr_1.2.1.tar.gz \
    MultiAssayExperiment_1.24.0.tar.gz \
    CAGEr_2.4.0.tar.gz
