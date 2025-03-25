FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies first
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

# Verify ps command is available
RUN which ps && ps --version

# Install required R packages
RUN R -e "install.packages(c('optparse', 'data.table', 'GenomicRanges', 'rtracklayer'), repos='https://cloud.r-project.org/')"

# Install SPLICE-Q
RUN pip3 install splice-q

# Set environment path to include system binaries
ENV PATH="/usr/sbin:/usr/bin:/sbin:/bin:${PATH}"

# Create a non-root user
RUN useradd -m -s /bin/bash spliceq && \
    usermod -aG staff spliceq

# Set working directory
WORKDIR /data

# Switch to non-root user
USER spliceq

# Test ps command as non-root user
RUN ps --version
