#!/bin/bash

# RNA-seq Pipeline Launcher with Automatic Docker Management
# This script ensures Docker is running before launching the Nextflow pipeline
# Usage: bash run_pipeline.sh [nextflow options]
# Example: bash run_pipeline.sh --single_end true --outdir my_results/

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if Docker daemon is running
check_docker_running() {
    if docker info >/dev/null 2>&1; then
        return 0  # Docker is running
    else
        return 1  # Docker is not running
    fi
}

# Function to start Docker
start_docker() {
    print_info "Docker is not running. Attempting to start Docker..."
    
    # Check if systemd is available (most modern Linux systems)
    if command -v systemctl >/dev/null 2>&1; then
        print_info "Using systemctl to start Docker..."
        
        # Try to start without sudo first (if user has permissions)
        if systemctl start docker >/dev/null 2>&1; then
            print_success "Docker started successfully"
            return 0
        fi
        
        # If that fails, try with sudo
        print_warning "Requires sudo privileges to start Docker"
        if sudo systemctl start docker; then
            print_success "Docker started successfully with sudo"
            return 0
        else
            print_error "Failed to start Docker with systemctl"
            return 1
        fi
        
    # Check for service command (older init systems)
    elif command -v service >/dev/null 2>&1; then
        print_info "Using service command to start Docker..."
        
        if sudo service docker start; then
            print_success "Docker started successfully"
            return 0
        else
            print_error "Failed to start Docker with service command"
            return 1
        fi
    
    else
        print_error "Could not find systemctl or service command"
        print_error "Please start Docker manually before running the pipeline"
        return 1
    fi
}

# Function to wait for Docker to be ready
wait_for_docker() {
    print_info "Waiting for Docker daemon to be ready..."
    local max_attempts=30
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        if check_docker_running; then
            print_success "Docker is ready"
            return 0
        fi
        
        echo -n "."
        sleep 1
        attempt=$((attempt + 1))
    done
    
    echo ""
    print_error "Docker daemon did not become ready in time"
    return 1
}

# Function to check Docker installation
check_docker_installed() {
    if ! command -v docker >/dev/null 2>&1; then
        print_error "Docker is not installed"
        print_info "Please run one of the setup scripts first:"
        print_info "  - bash prepare.sh (for Red Hat/CentOS/Fedora)"
        print_info "  - bash install_for_ami2023linux_aws.sh (for AWS AMI 2023 Linux)"
        exit 1
    fi
}

# Function to check Nextflow installation
check_nextflow_installed() {
    if ! command -v nextflow >/dev/null 2>&1; then
        print_error "Nextflow is not installed or not in PATH"
        print_info "Please run one of the setup scripts first, or add Nextflow to your PATH:"
        print_info "  export PATH=\"\$HOME/.local/bin:\$PATH\""
        exit 1
    fi
}

# Main execution starts here
print_info "RNA-seq Pipeline Launcher"
echo "========================================"

# Check if Docker is installed
check_docker_installed

# Check if Nextflow is installed
check_nextflow_installed

# Check if Docker is running
if check_docker_running; then
    print_success "Docker is already running"
else
    # Try to start Docker
    if start_docker; then
        # Wait for Docker to be fully ready
        if ! wait_for_docker; then
            exit 1
        fi
    else
        print_error "Could not start Docker automatically"
        print_info "Please start Docker manually with one of these commands:"
        print_info "  sudo systemctl start docker"
        print_info "  sudo service docker start"
        exit 1
    fi
fi

# Verify Docker is accessible
print_info "Verifying Docker access..."
if docker ps >/dev/null 2>&1; then
    print_success "Docker is accessible"
else
    print_error "Docker is running but current user cannot access it"
    print_info "You may need to:"
    print_info "  1. Add your user to the docker group: sudo usermod -aG docker \$USER"
    print_info "  2. Log out and log back in for the group change to take effect"
    print_info "  3. Or run this script with sudo"
    exit 1
fi

# Launch the Nextflow pipeline
echo ""
print_info "Launching RNA-seq pipeline with Nextflow..."
print_info "Pipeline arguments: $*"
echo "========================================"
echo ""

# Run the pipeline with all passed arguments
nextflow run main.nf "$@"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    print_success "Pipeline completed successfully!"
else
    echo ""
    print_error "Pipeline failed. Check the error messages above."
    exit 1
fi
