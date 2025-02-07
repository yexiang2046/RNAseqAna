FROM ubuntu:20.04

# Set environment variable to suppress interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    bzip2 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    bedtools \
    unzip \
    libisal-dev \
    libdeflate-dev \
    openjdk-8-jre \
    zlib1g \
    libstdc++6 \
    libbz2-1.0 \
    liblzma5 \
    libncurses5 \
    python3-pip

# Download and build htslib
RUN cd /usr/local/src && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && make \
    && make install

# Download and build samtools
RUN cd /usr/local/src && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -vxjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make \
    && make install

# Install MultiQC via pip
RUN pip3 install multiqc

# Install STAR
RUN wget -O /usr/local/bin/STAR.zip https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.zip \
    && unzip /usr/local/bin/STAR.zip -d /usr/local/bin/ \
    && cd /usr/local/bin/STAR-2.7.9a/source \
    && make STAR \
    && ln -s /usr/local/bin/STAR-2.7.9a/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR

# Install fastp
RUN wget -O /usr/local/bin/fastp.zip https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.2.zip \
    && unzip /usr/local/bin/fastp.zip -d /usr/local/bin/ \
    && cd /usr/local/bin/fastp-0.23.2 \
    && make \
    && ln -s /usr/local/bin/fastp-0.23.2/fastp /usr/local/bin/fastp

# Install subread
RUN wget -O /usr/local/bin/subread.tar.gz https://downloads.sourceforge.net/project/subread/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz \
    && tar -xzf /usr/local/bin/subread.tar.gz -C /usr/local/bin/ \
    && ln -s /usr/local/bin/subread-2.0.3-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts

# Install FastQC
RUN wget -O fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/FastQC/fastqc \
    && ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc

# Set environment variables
ENV PATH="/usr/local/htslib-1.9:${PATH}"
ENV PATH="/usr/local/samtools-1.9:${PATH}"

# Verify installations
RUN java -version && \
    fastqc --version && \
    multiqc --version && \
    STAR --version && \
    fastp --version && \
    featureCounts -v && \
    samtools --version

COPY . /usr/src/RNAseqAna
WORKDIR /usr/src/RNAseqAna