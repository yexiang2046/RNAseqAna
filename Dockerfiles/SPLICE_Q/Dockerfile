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

# Verify ps is installed and in PATH
RUN which ps

# Install required R packages
RUN R -e "install.packages(c('optparse', 'data.table', 'GenomicRanges', 'rtracklayer'), repos='https://cloud.r-project.org/')"

# Install SPLICE-Q
RUN pip3 install splice-q

# Create a non-root user and add to necessary groups
RUN useradd -m -s /bin/bash spliceq && \
    usermod -aG staff spliceq

# Set environment path to include system binaries
ENV PATH="/usr/bin:/bin:${PATH}"

# Switch to non-root user
USER spliceq

# Set working directory
WORKDIR /data

# Verify ps command is accessible for non-root user
RUN ps --version
