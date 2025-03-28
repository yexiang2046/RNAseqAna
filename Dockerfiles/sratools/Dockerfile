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
    && rm -rf /var/lib/apt/lists/*

# Clone and build SRA tools
WORKDIR /tmp
RUN git clone https://github.com/ncbi/sra-tools.git && \
    cd sra-tools && \
    git checkout 3.2.1 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf sra-tools

# Set up environment
ENV PATH="/usr/local/bin:${PATH}"

# Create a non-root user
RUN useradd -m -s /bin/bash sratools

# Switch to non-root user
USER sratools

# Set working directory
WORKDIR /home/sratools

# Default command
CMD ["/bin/bash"] 