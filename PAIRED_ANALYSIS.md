# Paired Sample Analysis Guide

## Overview
The RNAseq pipeline now supports differential expression analysis with paired samples using edgeR's paired design framework.

## When to Use Paired Analysis
Use paired analysis when samples are related, such as:
- Before/after treatment from same patient
- Matched tumor/normal tissue pairs
- Time series from same subject
- Technical replicates from same biological sample

## Metadata Format

### Unpaired Analysis (default)
```
SampleId    group
sample1     control
sample2     control
sample3     treated
sample4     treated
```

### Paired Analysis
Add a column identifying which samples are paired:
```
SampleId    group       patient_id
sample1     control     patient1
sample2     treated     patient1
sample3     control     patient2
sample4     treated     patient2
```

The paired column name can be anything (e.g., `patient_id`, `subject`, `pair_id`).

## Usage

### Running via Command Line
```bash
edger.r -c counts.txt -m metadata.txt -g annotation.gtf -o results -p patient_id
```

The `-p` flag specifies the column name for pairing.

### Running via Nextflow
Currently, the Nextflow module needs to be updated to support the paired parameter.

## Statistical Model

### Unpaired Design
- Formula: `~0 + group`
- Tests all pairwise group comparisons
- Assumes samples are independent

### Paired Design
- Formula: `~pair + group`
- Tests group effects while accounting for paired structure
- Blocks on the pairing variable to reduce within-pair variation
- More powerful when pairing is strong

## Output
- Same output files as unpaired analysis
- DEG_*.csv files contain differential expression results
- normalized_counts_CPM.csv contains normalized expression values
- PCA_plot.png visualizes sample relationships

## Example Scenarios

### Scenario 1: Pre/Post Treatment
```
SampleId    group       subject
pre_1       pre         patient1
post_1      post        patient1
pre_2       pre         patient2
post_2      post        patient2
```
Command: `edger.r ... -p subject`

### Scenario 2: Tumor-Normal Pairs
```
SampleId    group       patient
normal_1    normal      P001
tumor_1     tumor       P001
normal_2    normal      P002
tumor_2     tumor       P002
```
Command: `edger.r ... -p patient`

### Scenario 3: Multiple Time Points (use blocking)
```
SampleId    group       subject
t0_1        t0          S1
t24_1       t24         S1
t48_1       t48         S1
t0_2        t0          S2
t24_2       t24         S2
t48_2       t48         S2
```
Command: `edger.r ... -p subject`

## Notes
- Pairing variable must be balanced across groups
- Each pair should have samples in both conditions being compared
- Paired analysis requires at least 3 pairs for reliable results
- Paired design increases power when biological variation between pairs is high
