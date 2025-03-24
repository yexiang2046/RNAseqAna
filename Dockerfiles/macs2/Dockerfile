# Use continuumio/miniconda3 as base image
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Initialize conda and create environment
RUN conda init bash && \
    . ~/.bashrc && \
    conda create -n macs2_env python=3.9 -y && \
    echo "conda activate macs2_env" >> ~/.bashrc && \
    . ~/.bashrc && \
    conda install -c bioconda macs2=2.2.7.1 -y && \
    conda install -c conda-forge r-tidyverse r-ggplot2 r-dplyr r-tidyr -y

# Add conda environment to PATH
ENV PATH="/opt/conda/envs/macs2_env/bin:${PATH}"

# Create a non-root user
RUN useradd -m -u 1000 appuser && \
    chown -R appuser:appuser /app

# Switch to non-root user
USER appuser
