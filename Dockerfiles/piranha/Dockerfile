# Use Ubuntu as base image
FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    bedtools \
    wget \
    tar \
    unzip \
    libgsl-dev \
    python3 \
    python-is-python3 \
    python3-pip \
    git \
    cmake \
    pkg-config \
    zlib1g-dev \
    libjsoncpp-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

# install pandas, matplotlib, seaborn for generating reports
RUN pip3 install --no-cache-dir \
     pandas \
     matplotlib \
     seaborn

# Set working directory
WORKDIR /opt

# Install BAMTools first (required for Piranha BAM support)
RUN git clone https://github.com/pezmaster31/bamtools.git \
    && cd bamtools \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && ldconfig \
    && cd ../.. \
    && rm -rf bamtools

# Set environment variables for BAMTools
ENV BAMTOOLS_ROOT=/usr/local
ENV CPPFLAGS="-I${BAMTOOLS_ROOT}/include/bamtools"
ENV LDFLAGS="-L${BAMTOOLS_ROOT}/lib -Wl,-rpath,${BAMTOOLS_ROOT}/lib"

# Download and install Piranha with BAMTools support
RUN wget https://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz \
    && tar -xf piranha-1.2.1.tar.gz \
    && cd piranha-1.2.1 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf piranha-1.2.1.tar.gz

# Add Piranha to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf samtools-1.21.tar.bz2 \
    && cd samtools-1.21 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.21*

# Create symlink after installation
RUN if [ -f /usr/local/bin/piranha ]; then \
        ln -sf /usr/local/bin/piranha /usr/local/bin/Piranha; \
    fi

# Verify installations
RUN which piranha || echo "piranha not found" && \
    which samtools || echo "samtools not found"

# Set default working directory
WORKDIR /data