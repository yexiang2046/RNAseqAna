# Stage 1: Build stage
FROM ubuntu:20.04 as build

# Install necessary packages
RUN apt-get update && apt-get install -y build-essential wget bzip2 zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev bedtools

# Install ncurses development headers (required for htslib)
RUN apt-get install -y libncurses-dev

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

COPY . /usr/src/ASPEAK
WORKDIR /usr/src/ASPEAK

# Check if bedtools works (for debugging)
RUN bedtools --version || echo "Bedtools failed to run"