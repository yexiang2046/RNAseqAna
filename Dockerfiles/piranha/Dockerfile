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
    default-jre \
    default-jdk \
    gradle \
    && rm -rf /var/lib/apt/lists/*

# install pandas, matplotlib, seaborn for generating reports
RUN pip3 install --no-cache-dir \
    pandas \
    matplotlib \
    seaborn

# Set working directory
WORKDIR /opt

# Install htslib first (required for samtools)
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xjf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.21*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.20*

# Install Picard tools from source
RUN mkdir -p /opt/picard && \
    wget https://github.com/broadinstitute/picard/archive/refs/tags/3.3.0.tar.gz && \
    tar -xzf 3.3.0.tar.gz && \
    cd picard-3.3.0 && \
    ./gradlew shadowJar && \
    cp build/libs/picard.jar /opt/picard/ && \
    cd .. && \
    rm -rf picard-3.3.0 3.3.0.tar.gz

# Install BAMTools for BAM file support
RUN git clone https://github.com/pezmaster31/bamtools.git \
    && cd bamtools \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && ldconfig

# Download and install Piranha
RUN wget https://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz \
    && tar -xf piranha-1.2.1.tar.gz \
    && cd piranha-1.2.1 \
    && ./configure \
    && make all \
    && make install

# Add Piranha and Picard to PATH and set PICARD_HOME
ENV PATH="/opt/piranha-1.2.1/bin:${PATH}"
ENV PICARD_HOME="/opt/picard"