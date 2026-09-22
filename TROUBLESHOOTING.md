# Troubleshooting: Pipeline Termination After STAR_INDEX

## Quick Diagnostic

Run the diagnostic script from the directory where you executed `nextflow run main.nf`:

```bash
bash diagnose_pipeline.sh
```

## Common Issues and Solutions

### Issue 1: Missing Reference Genome File

**Symptom**: Pipeline stops immediately after starting, or STAR_INDEX fails

**Cause**: No `*.genome.fa` file found in the project directory

**Solution**:
```bash
# Download a reference genome (example: human GRCh38)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa hg38.genome.fa
```

### Issue 2: Insufficient Memory for STAR_INDEX

**Symptom**: STAR_INDEX fails with "out of memory" error, or container exits with error code 137

**Cause**: STAR genome indexing requires ~30-40GB RAM for human genome

**Solution**:
1. Ensure your system has at least 64GB RAM available
2. Check memory usage: `free -h`
3. Reduce memory for other processes
4. Consider using a pre-built STAR index:

```bash
# Use pre-built index
nextflow run main.nf --star_index /path/to/prebuilt/star_index/
```

### Issue 3: No Input FASTQ Files

**Symptom**: ALIGN process never starts after STAR_INDEX completes

**Cause**: No FASTQ files matching the expected patterns in `data/`

**For Paired-End (default)**:

The pipeline automatically detects common naming patterns:
- `sample_R1.fastq.gz` / `sample_R2.fastq.gz`
- `sample_1.fastq.gz` / `sample_2.fastq.gz`
- `sample.R1.fastq.gz` / `sample.R2.fastq.gz`
- `sampleR1.fastq.gz` / `sampleR2.fastq.gz`
- And more...

```bash
# Check your files:
ls data/*.fastq.gz
ls data/*.fq.gz

# If your files use a different pattern, specify it:
nextflow run main.nf --read_pattern '*_read{1,2}.fq.gz'
```

**For Single-End**:
```bash
# Run with single-end flag:
nextflow run main.nf --single_end true

# Files must match: *.fastq.gz or *.fq.gz
ls data/*.fastq.gz
```

### Issue 4: TRIM Process Failed

**Symptom**: STAR_INDEX completes but ALIGN never starts

**Cause**: fastp (TRIM process) failed to process input reads

**Solution**:
```bash
# Check TRIM process logs in work directory
find work -name ".command.sh" -exec grep -l "fastp" {} \; | head -1 | xargs dirname

# Check that directory's .command.err and .command.out files
```

### Issue 5: Docker Container Issues

**Symptom**: Process fails to start or exits immediately

**Cause**: Docker not running, or image pull failures

**Solution**:
```bash
# Check Docker is running
docker ps

# Manually pull required images
docker pull quay.io/biocontainers/star:2.7.11b--h5ca1c30_6
docker pull staphb/fastp:0.24.0
docker pull xiang2019/rnaseq_cmd:v1.0.0
docker pull multiqc/multiqc:pdf-v1.34
```

### Issue 6: Configuration Path Issues

**Symptom**: Pipeline can't find GTF or data files

**Cause**: Default paths in `main.nf` don't match your system

**Solution**:
Override parameters when running:
```bash
nextflow run main.nf \
  --data_dir /path/to/your/fastq/files \
  --gtf /path/to/your/annotation.gtf \
  --outdir my_results
```

## How to Read Nextflow Logs

### Main log file
```bash
# Full log
less .nextflow.log

# Filter for specific process
grep "STAR_INDEX" .nextflow.log
grep "ERROR" .nextflow.log
```

### Process-specific logs (in work directory)
Each task has a unique work directory under `work/XX/YYYYYYYY.../`

Files in each task directory:
- `.command.sh` - The script that was executed
- `.command.out` - Standard output
- `.command.err` - Standard error
- `.exitcode` - Exit status (0 = success)
- `.command.log` - Combined log

```bash
# Find a specific process's work directory
find work -name ".command.sh" -exec grep -l "STAR.*genomeGenerate" {} \; | head -1 | xargs dirname

# Then inspect that directory
cd <work_directory>
cat .exitcode      # Check if it succeeded (0) or failed
cat .command.err   # Check error messages
cat .command.out   # Check output
```

## Resume a Failed Run

After fixing the issue, resume from where it stopped:

```bash
nextflow run main.nf -resume
```

This will reuse successfully completed tasks and only re-run failed ones.

## Getting More Help

If the issue persists, collect this information:

1. Full Nextflow log: `cat .nextflow.log`
2. Failed process work directory contents
3. Output of diagnostic script
4. System specs: `free -h`, `nproc`, `docker --version`
5. Nextflow version: `nextflow -version`
