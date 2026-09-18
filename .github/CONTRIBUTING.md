# Contributing to RNAseqAna

Thank you for your interest in contributing to RNAseqAna! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Testing](#testing)
- [Coding Standards](#coding-standards)
- [Submitting Changes](#submitting-changes)
- [CI/CD Pipeline](#cicd-pipeline)

## Code of Conduct

This project adheres to a code of conduct that all contributors are expected to follow. Please be respectful and professional in all interactions.

## Getting Started

### Prerequisites

- Nextflow (>=21.04.0)
- Docker or compatible container runtime
- Git
- Basic knowledge of RNA-seq analysis workflows

### Setting Up Development Environment

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/RNAseqAna.git
   cd RNAseqAna
   ```

3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/ORIGINAL_OWNER/RNAseqAna.git
   ```

4. Install prerequisites:
   ```bash
   bash prepare.sh
   ```

## Development Workflow

### Creating a Feature Branch

1. Ensure your master branch is up to date:
   ```bash
   git checkout master
   git pull upstream master
   ```

2. Create a new feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

### Making Changes

1. **Nextflow Pipeline Changes**: Edit files in `main.nf` or `modules/`
2. **R Script Changes**: Edit files in `bin/` or `R_scripts/`
3. **Documentation**: Update `README.md`, `CLAUDE.md`, or other docs as needed

### Commit Guidelines

Use clear, descriptive commit messages:

```
<type>: <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks
- `ci`: CI/CD changes

**Example:**
```
feat: add support for single-end reads in STAR alignment

- Updated ALIGN process to handle single-end mode
- Added conditional logic for read file inputs
- Updated documentation with single-end examples

Closes #123
```

## Testing

### Running Tests Locally

#### Nextflow Pipeline Tests

```bash
# Validate configuration
nextflow config -profile standard

# Run with test data (requires STAR index and GTF)
bash tests/run_tests.sh --star_index /path/to/index --gtf /path/to/annotation.gtf
```

#### R Script Tests

```bash
# Test edgeR script
docker run --rm \
  -v $(pwd):/data \
  xiang2019/rnaseq_renv:v1.0.2 \
  Rscript /data/bin/edger.r --help

# Test with sample data
docker run --rm \
  -v $(pwd):/data \
  xiang2019/rnaseq_renv:v1.0.2 \
  Rscript /data/bin/edger.r \
    -c /data/test_counts.txt \
    -m /data/test_metadata.txt \
    -g /data/test_annotation.gtf \
    -o /data/test_output \
    -s human
```

### Writing Tests

- Add test cases for new features
- Ensure tests are reproducible
- Use synthetic or public data for tests
- Document test data requirements

## Coding Standards

### Nextflow

- Use descriptive process names in UPPERCASE
- Include comments for complex logic
- Pin container versions (avoid `:latest`)
- Set appropriate resource directives (cpus, memory, time)
- Use `params` for configurable values
- Follow Nextflow DSL2 syntax

**Example:**
```groovy
process ALIGN {
    container 'quay.io/biocontainers/star:2.7.11b'
    cpus 16
    memory '64 GB'
    time '12h'

    input:
    path star_index
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    STAR --genomeDir ${star_index} \
         --readFilesIn ${reads} \
         --outSAMtype BAM SortedByCoordinate
    """
}
```

### R Scripts

- Use clear variable names
- Include help/usage information
- Handle errors gracefully
- Use existing Bioconductor packages when possible
- Comment complex statistical procedures
- Follow tidyverse style guide where applicable

**Example:**
```r
#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type="character",
              help="Path to count matrix file"),
  make_option(c("-m", "--metadata"), type="character",
              help="Path to metadata file")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Validate inputs
if (is.null(opt$counts)) {
  stop("Count matrix file is required")
}
```

### Documentation

- Keep README.md up to date
- Document all parameters and options
- Include usage examples
- Update CLAUDE.md for major architectural changes
- Add comments for non-obvious code

## Submitting Changes

### Before Submitting

1. **Test your changes:**
   ```bash
   # Run local tests
   nextflow config -profile standard
   
   # Run R script validation
   docker run --rm -v $(pwd):/data xiang2019/rnaseq_renv:v1.0.2 \
     Rscript -e "source('/data/bin/your_script.r')"
   ```

2. **Update documentation:**
   - README.md for user-facing changes
   - CLAUDE.md for architectural changes
   - Code comments for complex logic

3. **Commit your changes:**
   ```bash
   git add .
   git commit -m "feat: your descriptive message"
   ```

4. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

### Creating a Pull Request

1. Go to your fork on GitHub
2. Click "New Pull Request"
3. Select your feature branch
4. Fill out the PR template completely
5. Ensure CI checks pass

### PR Review Process

1. **Automated Checks**: All CI/CD workflows must pass
2. **Code Review**: At least one maintainer will review your code
3. **Testing**: Changes will be tested with sample data
4. **Documentation Review**: Ensure documentation is complete
5. **Approval**: Maintainer approval required before merge

## CI/CD Pipeline

This project uses GitHub Actions for automated testing:

### Workflows

- **CI** (`ci.yml`): Runs on every push/PR
  - Lints Nextflow code
  - Validates R scripts
  - Checks container availability
  - Tests pipeline structure

- **PR Checks** (`pr-checks.yml`): Additional PR validation
  - PR title validation
  - Changed file analysis
  - Size labeling
  - Breaking change detection

- **Test R Scripts** (`test-r-scripts.yml`): R-specific tests
  - Syntax validation
  - Dependency checks
  - Script testing with sample data

- **Security** (`security.yml`): Security scanning
  - Container vulnerability scanning
  - Secret detection
  - Dependency review
  - Dockerfile linting

- **Nightly** (`nightly.yml`): Scheduled comprehensive tests
  - Extended compatibility tests
  - Container update checks
  - Performance baselines

- **Release** (`release.yml`): Release automation
  - Version validation
  - Changelog generation
  - GitHub release creation

### Local CI Testing

You can test many CI checks locally before pushing:

```bash
# Validate Nextflow
nextflow config -profile standard

# Test R script syntax
docker run --rm -v $(pwd):/data xiang2019/rnaseq_renv:v1.0.2 \
  Rscript -e "parse('bin/your_script.r')"

# Check for secrets (install gitleaks)
gitleaks detect --source . --verbose

# Lint Dockerfiles (install hadolint)
hadolint Dockerfiles/*
```

## Questions or Problems?

- Open an issue for bugs or feature requests
- Use discussions for questions
- Tag maintainers for urgent issues

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (see LICENSE.md).

---

Thank you for contributing to RNAseqAna! 🧬
