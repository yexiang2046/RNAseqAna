#!/bin/bash

# This script prepares the environment for running the RNAseqAna pipeline on Red Hat-based systems (e.g., Fedora, CentOS, RHEL).
# It assumes a modern Red Hat system using dnf (Fedora/RHEL 8+). For older yum-based systems, replace 'dnf' with 'yum' where applicable.
# Requires sudo privileges for package installations.
# Run this script with: bash install.sh

set -e  # Exit on error

echo "Updating package lists..."
sudo dnf update -y

echo "Installing prerequisites for SDKMAN! (curl, zip, unzip)..."
sudo dnf install -y curl zip unzip

echo "Installing SDKMAN! for Java management..."
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"

echo "Installing Java via SDKMAN!..."
sdk install java 17.0.10-tem

echo "Installing Git (if not already installed)..."
sudo dnf install -y git

echo "Cloning the RNAseqAna repository..."
git clone https://github.com/yexiang2046/RNAseqAna.git
cd RNAseqAna

echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
export PATH="$HOME/.local/bin:$PATH"

echo "Installing Docker..."
sudo dnf install -y docker
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker $USER  # Add current user to docker group (log out and back in for changes to take effect)

echo "Installing R..."
sudo dnf install -y R

echo "Installing required R packages..."
Rscript -e '
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
}
BiocManager::install(c("edgeR", "limma"), ask=FALSE)
install.packages(c("pheatmap", "RColorBrewer", "optparse", "tidyverse"), repos="http://cran.us.r-project.org")
'

echo "Setup complete! You may need to log out and log back in for Docker group changes and to ensure PATH updates (for Nextflow)."
echo "To run the pipeline, use: nextflow run main.nf [options]"
