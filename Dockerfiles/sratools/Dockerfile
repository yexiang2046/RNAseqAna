FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
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

# Install SRA toolkit using the official binaries
WORKDIR /tmp
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    cp -r sratoolkit.*/bin/* /usr/local/bin/ && \
    rm -rf sratoolkit.*

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