FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    flex \
    bison \
    default-jdk \
    python3 \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Build ncbi-vdb and sra-tools
WORKDIR /tmp
RUN git clone https://github.com/ncbi/ncbi-vdb.git && \
    git clone https://github.com/ncbi/sra-tools.git && \
    cd ncbi-vdb && \
    ./configure && \
    make && \
    make install && \
    cd ../sra-tools && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf ncbi-vdb sra-tools

# Set up environment
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Create a non-root user
RUN useradd -m -s /bin/bash sratools

# Switch to non-root user
USER sratools

# Set working directory
WORKDIR /home/sratools

# Default command
CMD ["/bin/bash"] 