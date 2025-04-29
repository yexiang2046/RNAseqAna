# Use debian:bullseye-slim as base image
FROM debian:bullseye-slim

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
    wget \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba -O > /usr/local/bin/micromamba && \
    chmod u+x /usr/local/bin/micromamba

# Create micromamba environment
RUN micromamba create -n macs2_env -c bioconda -c conda-forge \
    python=3.9 \
    macs2=2.2.7.1 \
    r-tidyverse \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    -y

# Initialize micromamba shell
RUN micromamba shell init --shell=bash && \
    echo "micromamba activate macs2_env" >> ~/.bashrc

# Create entrypoint script
RUN echo '#!/bin/bash\nsource ~/.bashrc\nexec "$@"' > /entrypoint.sh && \
    chmod +x /entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]
CMD ["macs2"]



