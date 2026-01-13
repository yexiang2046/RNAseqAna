#!/bin/bash
# Build script for rnaseq_renv Docker container

set -e

VERSION="v1.0.2"
IMAGE_NAME="xiang2019/rnaseq_renv"

echo "Building ${IMAGE_NAME}:${VERSION}..."

# Copy the latest R scripts from bin directory
echo "Copying latest R scripts..."
mkdir -p bin
cp ../../bin/edger.r bin/
cp ../../bin/functional_analysis.r bin/

# Build the Docker image
echo "Building Docker image..."
docker build -t ${IMAGE_NAME}:${VERSION} .

# Tag as latest
echo "Tagging as latest..."
docker tag ${IMAGE_NAME}:${VERSION} ${IMAGE_NAME}:latest

echo ""
echo "Build complete!"
echo ""
echo "To push to Docker Hub, run:"
echo "  docker push ${IMAGE_NAME}:${VERSION}"
echo "  docker push ${IMAGE_NAME}:latest"
echo ""
echo "To update nextflow.config, change:"
echo "  container = '${IMAGE_NAME}:v1.0.0' → container = '${IMAGE_NAME}:${VERSION}'"
