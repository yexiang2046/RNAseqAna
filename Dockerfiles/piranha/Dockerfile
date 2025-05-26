# Use Ubuntu as base image (64-bit required by Piranha)
FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
# - GCC > 4.1 required for TR1 support
# - GSL required for scientific computing
# - Other build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    make \
    wget \
    tar \
    git \
    cmake \
    libgsl-dev \
    zlib1g-dev \
    autoconf \
    automake \
    libtool \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Clone and install BAMTools (required for BAM support)
RUN git clone https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    # Create symlinks for headers
    ln -s /usr/local/include/bamtools/api /usr/local/include/api && \
    cd ../.. && \
    rm -rf bamtools

# Set environment variables for BAMTools
ENV BAMTOOLS_ROOT=/usr/local
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV CPATH=/usr/local/include/bamtools:$CPATH

# Download and install Piranha
RUN wget https://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz && \
    tar -xf piranha-1.2.1.tar.gz && \
    cd piranha-1.2.1 && \
    # Initialize build system
    autoreconf -i && \
    # Configure with BAMTools support
    ./configure \
        CPPFLAGS="-I/usr/local/include/bamtools -I/usr/local/include" \
        LDFLAGS="-L/usr/local/lib" \
        LIBS="-lbamtools -lz" && \
    # Build and install
    make CPPFLAGS="-I/usr/local/include/bamtools -I/usr/local/include" && \
    make install && \
    # Cleanup
    cd .. && \
    rm -rf piranha-1.2.1*

# Verify installation and BAM support
RUN Piranha --help | grep -i "bam" || echo "BAM support not found"

# Set working directory for data
WORKDIR /data