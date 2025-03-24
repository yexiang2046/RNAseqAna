# Use Python 3.9 slim image as base
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    macs2==2.2.7.1 \
    numpy \
    pandas \
    scipy

# Install R packages
RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'dplyr', 'tidyr'), repos='https://cloud.r-project.org/')"

# Set environment variables
ENV PATH="/app:${PATH}"

# Create a non-root user
RUN useradd -m -u 1000 appuser && \
    chown -R appuser:appuser /app

# Switch to non-root user
USER appuser
