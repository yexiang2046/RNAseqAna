FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install gtftools
RUN pip install --no-cache-dir argparse sys numpy gtftools

# Set working directory
WORKDIR /data

# Add version and usage information
LABEL maintainer="Xiang Ye <yexiang2046@gmail.com>"
LABEL version="1.0"
LABEL description="Docker image for gtftools"