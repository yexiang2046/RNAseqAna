# Stage 1: Build stage
FROM ubuntu:20.04 AS build

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bedtools \
    libbz2-1.0 \
    libbz2-dev \
    liblzma5 \
    liblzma-dev \
    libncursesw6 \
    && rm -rf /var/lib/apt/lists/*

# Install samtools and htslib
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Build htslib
RUN cd /usr/local/src \
    && wget https://github.com/samtools/htslib/archive/refs/tags/1.9.tar.gz \
    && tar -xzf 1.9.tar.gz \
    && rm 1.9.tar.gz \
    && cd htslib-1.9 \
    && make \
    && make install

# Build samtools
RUN cd /usr/local/src \
    && wget https://github.com/samtools/samtools/archive/refs/tags/1.9.tar.gz \
    && tar -xzf 1.9.tar.gz \
    && rm 1.9.tar.gz \
    && cd samtools-1.9 \
    && ln -s /usr/local/src/htslib-1.9 htslib \
    && make \
    && make install

# Stage 2: Final stage
FROM quay.io/biocontainers/perl-math-cdf:0.1--pl5321h7b50bb2_11

# Copy the necessary libraries from the build stage
COPY --from=build /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libz.so.1
COPY --from=build /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6

# Copy the built binaries from the build stage
COPY --from=build /usr/local/src/htslib-1.9 /usr/local/htslib-1.9
COPY --from=build /usr/local/src/samtools-1.9 /usr/local/samtools-1.9
COPY --from=build /usr/bin/bedtools /usr/bin/bedtools

# Copy necessary libraries for bedtools
COPY --from=build /lib/x86_64-linux-gnu/libbz2.so.1 /lib/x86_64-linux-gnu/libbz2.so.1
COPY --from=build /lib/x86_64-linux-gnu/liblzma.so.5 /lib/x86_64-linux-gnu/liblzma.so.5
COPY --from=build /lib/x86_64-linux-gnu/libncursesw.so.6 /lib/x86_64-linux-gnu/libncursesw.so.6

# Add the binaries to the PATH
ENV PATH="/usr/local/htslib-1.9:${PATH}"
ENV PATH="/usr/local/samtools-1.9:${PATH}"
ENV PATH="/usr/bin:${PATH}"

# Copy the ASPeak source code
COPY . /usr/src/ASPEAK
WORKDIR /usr/src/ASPEAK

# Check if bedtools works (for debugging)
RUN bedtools --version || echo "Bedtools failed to run"