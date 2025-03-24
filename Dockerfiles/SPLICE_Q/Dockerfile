FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    procps \
    python3 \
    python3-pip \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('optparse', 'data.table', 'GenomicRanges', 'rtracklayer'), repos='https://cloud.r-project.org/')"

# Install SPLICE-Q
RUN pip3 install splice-q

# Create a non-root user
RUN useradd -m -s /bin/bash spliceq
USER spliceq

# Set working directory
WORKDIR /data