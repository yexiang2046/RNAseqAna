---
name: Bug Report
about: Report a bug in the RNA-seq pipeline
title: '[BUG] '
labels: bug
assignees: ''
---

## Bug Description

<!-- A clear and concise description of what the bug is -->

## To Reproduce

Steps to reproduce the behavior:

1. Run command: `...`
2. With parameters: `...`
3. Using data: `...`
4. See error

## Expected Behavior

<!-- A clear and concise description of what you expected to happen -->

## Actual Behavior

<!-- What actually happened -->

## Error Messages

```
<!-- Paste error messages or logs here -->
```

## Environment

**Nextflow:**
- Version: 
- Run command: 

**Docker:**
- Version: 
- Images used: 

**System:**
- OS: 
- Available RAM: 
- Available CPUs: 

**Pipeline:**
- Branch/commit: 
- Profile: 
- Modified config: (yes/no)

## Input Data

- Data type: (paired-end/single-end)
- Number of samples: 
- FASTQ size (per file): 
- Genome reference: 

## Workflow Stage

<!-- Which stage did the pipeline fail at? -->

- [ ] TRIM (fastp)
- [ ] STAR_INDEX
- [ ] ALIGN
- [ ] FEATURECOUNT
- [ ] MULTIQC
- [ ] edgeR analysis
- [ ] Functional analysis
- [ ] Heatmap generation
- [ ] Other (specify)

## Additional Context

<!-- Add any other context about the problem here -->

## Possible Solution

<!-- Optional: suggest a fix or reason for the bug -->

## Related Issues

<!-- Link any related issues -->
