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

# Build ncbi-vdb first
WORKDIR /tmp
RUN git clone https://github.com/ncbi/ncbi-vdb.git && \
    cd ncbi-vdb && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_SHARED_LIBS=ON \
          .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf ncbi-vdb

# Clone and build SRA tools
RUN git clone https://github.com/ncbi/sra-tools.git && \
    cd sra-tools && \
    git checkout 3.2.1 && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_TOOLS_INTERNAL=OFF \
          -DBUILD_TOOLS_LOADERS=OFF \
          -DBUILD_TOOLS_TEST_TOOLS=OFF \
          -DTOOLS_ONLY=ON \
          -DVERSION=3.2.1 \
          -DVDB_INCDIR=/usr/local/include \
          -DVDB_LIBDIR=/usr/local/lib \
          .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf sra-tools

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