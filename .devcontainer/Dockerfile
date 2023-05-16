FROM bioconductor/bioconductor_docker:RELEASE_3_17

# install R packages needed for VSCode interaction and package management
RUN R -e "install.packages(c('languageserver', 'renv', 'remotes', 'httpgd'))"

# install R packages for data science
RUN R -e 'install.packages(c("devtools", "viridis", "tidyverse", "hdf5r"))'


# install Bioinformatic tools
RUN R -e 'install.packages(c("pheatmap", "enrichR", "BiocManager"))'
RUN R -e 'BiocManager::install(c("Rsubread", "edgeR", "limma", "msigdbr"))'



# install radian via python and pip3
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends python3-setuptools \
    rna-star \
    htop

RUN pip3 install radian

