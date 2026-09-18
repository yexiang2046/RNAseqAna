# CI/CD Pipeline Documentation

This document describes the Continuous Integration and Continuous Deployment (CI/CD) pipeline for the RNAseqAna project.

## Overview

The pipeline uses GitHub Actions to automate testing, validation, security scanning, and release management. All workflows are defined in `.github/workflows/`.

## Workflows

### 1. CI (`ci.yml`)

**Triggers:** Push to master/main, Pull Requests
**Purpose:** Core validation and testing

#### Jobs:

- **lint**: Validates Nextflow pipeline syntax and configuration
- **container-check**: Verifies all required Docker containers are available
- **validate-r-scripts**: Checks R script syntax in containers
- **test-pipeline**: Tests pipeline with synthetic data
- **documentation**: Validates README and documentation
- **summary**: Aggregates results and reports status

**Badge:**
```markdown
![CI](https://github.com/yourusername/RNAseqAna/workflows/CI/badge.svg)
```

### 2. PR Checks (`pr-checks.yml`)

**Triggers:** Pull Request events (opened, synchronized, reopened, ready_for_review)
**Purpose:** Additional PR validation and labeling

#### Jobs:

- **pr-validation**: Validates PR title and checks for breaking changes
- **changed-files**: Analyzes changed files and suggests reviewers
- **size-label**: Automatically labels PR size (XS, S, M, L, XL)
- **test-on-pr**: Tests configuration loading

**Features:**
- Automatic size labeling based on lines changed
- Breaking change detection
- Smart reviewer suggestions based on file types

### 3. Test R Scripts (`test-r-scripts.yml`)

**Triggers:** Changes to R scripts, manual dispatch
**Purpose:** Comprehensive R script testing

#### Jobs:

- **lint-r-scripts**: Syntax validation for all R scripts
- **test-edger**: Tests edgeR differential expression script
- **test-functional-analysis**: Tests functional enrichment script
- **test-heatmap**: Tests heatmap generation script
- **verify-containers**: Verifies R container dependencies

**Badge:**
```markdown
![Test R Scripts](https://github.com/yourusername/RNAseqAna/workflows/Test%20R%20Scripts/badge.svg)
```

### 4. Security (`security.yml`)

**Triggers:** Push, PR, Weekly schedule (Monday 00:00 UTC), Manual dispatch
**Purpose:** Security scanning and vulnerability detection

#### Jobs:

- **container-scan**: Scans configurations with Trivy
- **secret-scan**: Detects secrets with Gitleaks
- **dependency-review**: Reviews dependency changes in PRs
- **dockerfile-lint**: Lints Dockerfiles with hadolint
- **verify-container-sources**: Validates container image sources
- **permissions-check**: Checks file permissions and sensitive files
- **code-quality**: Scans for hardcoded credentials

**Badge:**
```markdown
![Security](https://github.com/yourusername/RNAseqAna/workflows/Security%20Scan/badge.svg)
```

### 5. Nightly Tests (`nightly.yml`)

**Triggers:** Daily at 2:00 AM UTC, Manual dispatch
**Purpose:** Extended testing and maintenance checks

#### Jobs:

- **extended-tests**: Tests with multiple Nextflow versions (23.10.0, 24.04.0, latest)
- **container-updates**: Checks for outdated container images
- **test-documentation**: Validates documentation examples
- **dependency-health**: Checks dependency versions and health
- **performance-baseline**: Measures pipeline startup time
- **notify-results**: Aggregates and reports nightly test results

**Badge:**
```markdown
![Nightly](https://github.com/yourusername/RNAseqAna/workflows/Nightly%20Tests/badge.svg)
```

### 6. Documentation (`docs.yml`)

**Triggers:** Changes to markdown files, Manual dispatch
**Purpose:** Documentation quality checks

#### Jobs:

- **lint-markdown**: Lints markdown files with markdownlint
- **check-links**: Validates links in documentation
- **validate-examples**: Checks code examples for safety
- **check-structure**: Verifies required documentation exists
- **spelling**: Spell-checking with misspell

### 7. Release (`release.yml`)

**Triggers:** Version tags (v*.*.*), Manual dispatch
**Purpose:** Automated release management

#### Jobs:

- **create-release**: Creates GitHub release with changelog
- **validate-release**: Validates release artifacts
- **notify**: Reports release status

**Badge:**
```markdown
![Release](https://github.com/yourusername/RNAseqAna/workflows/Release/badge.svg)
```

## Status Badges

Add these badges to your README.md:

```markdown
![CI](https://github.com/yourusername/RNAseqAna/workflows/CI/badge.svg)
![Security](https://github.com/yourusername/RNAseqAna/workflows/Security%20Scan/badge.svg)
![Test R Scripts](https://github.com/yourusername/RNAseqAna/workflows/Test%20R%20Scripts/badge.svg)
![Nightly](https://github.com/yourusername/RNAseqAna/workflows/Nightly%20Tests/badge.svg)
![Documentation](https://github.com/yourusername/RNAseqAna/workflows/Documentation/badge.svg)
```

## Local Testing

### Before Pushing

Run these checks locally to catch issues early:

```bash
# 1. Validate Nextflow config
nextflow config -profile standard

# 2. Test R script syntax
docker run --rm -v $(pwd):/data xiang2019/rnaseq_renv:v1.0.2 \
  Rscript -e "source('bin/edger.r', echo=FALSE)"

# 3. Check for secrets (requires gitleaks)
gitleaks detect --source . --verbose

# 4. Lint Dockerfiles (requires hadolint)
hadolint Dockerfiles/*

# 5. Check markdown (requires markdownlint-cli)
markdownlint '**/*.md' --ignore node_modules
```

### Install Local Tools

```bash
# Install gitleaks
brew install gitleaks  # macOS
# or download from https://github.com/gitleaks/gitleaks/releases

# Install hadolint
brew install hadolint  # macOS
# or download from https://github.com/hadolint/hadolint/releases

# Install markdownlint
npm install -g markdownlint-cli
```

## Workflow Configuration

### Secrets Required

No secrets are required for basic CI/CD. Optional secrets for enhanced features:

- `GITHUB_TOKEN`: Automatically provided by GitHub Actions
- Custom registry credentials (if using private containers)

### Branch Protection Rules

Recommended settings for the `master` branch:

- ✅ Require pull request reviews before merging
- ✅ Require status checks to pass before merging
  - Required checks: `lint`, `test-pipeline`, `validate-r-scripts`
- ✅ Require branches to be up to date before merging
- ✅ Include administrators
- ✅ Restrict who can push to matching branches

### PR Labels

Automatically applied labels:

- `size/XS`: < 10 lines changed
- `size/S`: 10-49 lines changed
- `size/M`: 50-249 lines changed
- `size/L`: 250-999 lines changed
- `size/XL`: 1000+ lines changed

Manual labels for organization:

- `bug`: Bug fix
- `enhancement`: New feature
- `documentation`: Documentation update
- `dependencies`: Dependency updates
- `breaking-change`: Breaking change
- `help-wanted`: Help wanted
- `good-first-issue`: Good for newcomers

## Dependabot

Automated dependency updates configured in `.github/dependabot.yml`:

- **GitHub Actions**: Weekly updates on Monday
- **Docker**: Weekly container base image checks

## Troubleshooting

### Workflow Fails

1. **Check the logs**: Click on the failed job in the Actions tab
2. **Run locally**: Use commands from "Local Testing" section
3. **Check requirements**: Ensure all required files exist

### Container Pull Failures

```bash
# Manually test container pull
docker pull xiang2019/rnaseq_cmd:v1.0.0
```

### R Script Failures

```bash
# Test R script in container
docker run --rm -v $(pwd):/data xiang2019/rnaseq_renv:v1.0.2 \
  Rscript /data/bin/your_script.r --help
```

### Nextflow Config Issues

```bash
# Validate config
nextflow config -profile standard -show-profiles

# Check specific process
nextflow config -profile standard | grep -A 10 "withName: ALIGN"
```

## Performance

### Workflow Timing (Approximate)

- **CI**: 5-8 minutes
- **PR Checks**: 2-3 minutes
- **Test R Scripts**: 8-12 minutes
- **Security**: 6-10 minutes
- **Documentation**: 2-4 minutes
- **Nightly**: 15-25 minutes

### Optimization Tips

1. **Use matrix builds** for parallel testing
2. **Cache dependencies** where possible
3. **Skip unnecessary jobs** with path filters
4. **Use self-hosted runners** for faster builds

## Monitoring

### GitHub Actions Usage

Monitor workflow usage in repository Settings → Actions → General

- View workflow run history
- Check concurrent workflow limits
- Monitor storage usage

### Notifications

Configure notifications in your GitHub settings:

- Email: Workflow failures
- Slack: Integration with GitHub Actions
- Status checks: Required for merges

## Best Practices

1. **Write descriptive commit messages** for automatic changelog
2. **Test locally** before pushing
3. **Keep workflows DRY** using reusable workflows
4. **Pin action versions** for stability
5. **Use matrix builds** for version testing
6. **Monitor workflow performance** regularly
7. **Update dependencies** via Dependabot PRs

## Future Enhancements

- [ ] Add performance benchmarking
- [ ] Implement code coverage reporting
- [ ] Add integration tests with real data
- [ ] Create deployment workflows for containers
- [ ] Add workflow visualization
- [ ] Implement artifact caching
- [ ] Add scheduled dependency updates

## References

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Nextflow CI/CD Best Practices](https://www.nextflow.io/docs/latest/sharing.html)
- [nf-core Tools](https://nf-co.re/tools)

---

For questions or suggestions about the CI/CD pipeline, please open an issue or discussion.
