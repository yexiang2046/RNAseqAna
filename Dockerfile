# Stage 1: Build stage
FROM ubuntu:20.04 as build

# Install necessary packages
RUN apt-get update && apt-get install -y build-essential wget bzip2 zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev bedtools unzip libisal-dev libdeflate-dev

# Install additional libraries
RUN apt-get install -y zlib1g libstdc++6

# Download and build htslib
RUN cd /usr/local/src && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && make

# Download and build samtools
RUN cd /usr/local/src && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -vxjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make

# Install FastQC
RUN wget -O /usr/local/bin/fastqc.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip /usr/local/bin/fastqc.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/FastQC/fastqc \
    && ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc

# Install MultiQC
RUN wget -O /usr/local/bin/multiqc.zip https://github.com/ewels/MultiQC/archive/refs/tags/v1.11.zip \
    && unzip /usr/local/bin/multiqc.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/MultiQC-1.11/multiqc \
    && ln -s /usr/local/bin/MultiQC-1.11/multiqc /usr/local/bin/multiqc

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

# Stage 2: Final stage
FROM quay.io/biocontainers/perl-math-cdf:0.1--pl5321h7b50bb2_11

# Copy the necessary libraries from the build stage
COPY --from=build /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libz.so.1
COPY --from=build /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6

# Copy the built binaries from the build stage
COPY --from=build /usr/local/src/htslib-1.9 /usr/local/htslib-1.9
COPY --from=build /usr/local/src/samtools-1.9 /usr/local/samtools-1.9
COPY --from=build /usr/bin/bedtools /usr/bin/bedtools
COPY --from=build /usr/local/bin/FastQC /usr/local/bin/FastQC
COPY --from=build /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc
COPY --from=build /usr/local/bin/MultiQC-1.11 /usr/local/bin/MultiQC-1.11
COPY --from=build /usr/local/bin/multiqc /usr/local/bin/multiqc
COPY --from=build /usr/local/bin/STAR-2.7.9a /usr/local/bin/STAR-2.7.9a
COPY --from=build /usr/local/bin/STAR /usr/local/bin/STAR
COPY --from=build /usr/local/bin/fastp-0.23.2 /usr/local/bin/fastp-0.23.2
COPY --from=build /usr/local/bin/fastp /usr/local/bin/fastp
COPY --from=build /usr/local/bin/subread-2.0.3-Linux-x86_64 /usr/local/bin/subread-2.0.3-Linux-x86_64
COPY --from=build /usr/local/bin/featureCounts /usr/local/bin/featureCounts

# Set environment variables
ENV PATH="/usr/local/htslib-1.9:${PATH}"
ENV PATH="/usr/local/samtools-1.9:${PATH}"

# Verify installations
RUN chmod +x /usr/local/bin/FastQC/fastqc && \
    chmod +x /usr/local/bin/fastqc && \
    chmod +x /usr/local/bin/multiqc && \
    chmod +x /usr/local/bin/STAR && \
    chmod +x /usr/local/bin/fastp && \
    chmod +x /usr/local/bin/featureCounts && \
    fastqc --version && \
    multiqc --version && \
    STAR --version && \
    fastp --version && \
    featureCounts -v && \
    samtools --version

COPY . /usr/src/RNAseqAna
WORKDIR /usr/src/RNAseqAna