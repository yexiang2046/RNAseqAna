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
     && rm -rf /var/lib/apt/lists/*
# install pandas, matplotlib, seaborn for generating reports
RUN pip3 install --no-cache-dir \
     pandas \
     matplotlib \
     seaborn
     

 # Set working directory
 WORKDIR /opt
 
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
 
 # Add Piranha to PATH
 ENV PATH="/opt/piranha-1.2.1/bin:${PATH}"