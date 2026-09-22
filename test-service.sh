#!/bin/bash

set -e

echo "=== RNA-seq Web Service Test Script ==="
echo ""

echo "1. Checking prerequisites..."
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed"
    exit 1
fi

if ! command -v docker-compose &> /dev/null; then
    echo "ERROR: Docker Compose is not installed"  
    exit 1
fi

echo "   ✓ Docker is installed"
echo "   ✓ Docker Compose is installed"
echo ""

echo "2. Validating Docker Compose configuration..."
if docker-compose config > /dev/null 2>&1; then
    echo "   ✓ Docker Compose configuration is valid"
else
    echo "   ✗ Docker Compose configuration has errors"
    exit 1
fi
echo ""

echo "=== All checks passed! ==="
echo ""
echo "To start the service, run:"
echo "  docker-compose up -d"
echo ""
