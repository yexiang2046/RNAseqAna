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
    autoconf \
    automake \
    libtool \
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
ENV CPPFLAGS="-I/usr/local/include/bamtools"
ENV LDFLAGS="-L/usr/local/lib -lbamtools"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Download and install Piranha with BAMTools support
RUN wget https://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz \
    && tar -xf piranha-1.2.1.tar.gz \
    && cd piranha-1.2.1 \
    && sed -i 's/AM_INIT_AUTOMAKE/AM_INIT_AUTOMAKE([subdir-objects])/' configure.ac \
    && autoreconf --force --install \
    && ./configure CPPFLAGS="-I/usr/local/include/bamtools -DWITH_BAM_SUPPORT" \
                  LDFLAGS="-L/usr/local/lib -lbamtools" \
                  LIBS="-lbamtools" \
    && make \
    && make install \
    && cd .. \
    && rm -rf piranha-1.2.1.tar.gz

# Add Piranha to PATH - use both possible installation locations
ENV PATH="/usr/local/bin:/opt/piranha-1.2.1/bin:${PATH}"

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf samtools-1.21.tar.bz2 \
    && cd samtools-1.21 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.21*

# Create symlink for Piranha
RUN ln -sf $(which piranha) /usr/local/bin/Piranha

# Verify installations and BAM support
RUN Piranha --help | grep -i "bam" || echo "BAM support not found"

# Set default working directory
WORKDIR /data