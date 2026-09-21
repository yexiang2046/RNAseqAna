#!/usr/bin/env python3
import gzip

# Create matching paired-end reads with exactly 75bp length
sequences_r1 = [
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
    "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
]

sequences_r2 = [
    "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
    "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG",
    "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
]

# Quality scores (same length as sequences)
quality = "I" * 75

# Write R1 file
with gzip.open("test_data/sample1_R1_001.fastq.gz", "wt") as f:
    for i, seq in enumerate(sequences_r1, 1):
        assert len(seq) == len(quality), f"Seq {i} length mismatch: {len(seq)} != {len(quality)}"
        f.write(f"@SEQ{i}\n")
        f.write(f"{seq}\n")
        f.write("+\n")
        f.write(f"{quality}\n")

# Write R2 file
with gzip.open("test_data/sample1_R2_001.fastq.gz", "wt") as f:
    for i, seq in enumerate(sequences_r2, 1):
        assert len(seq) == len(quality), f"Seq {i} length mismatch: {len(seq)} != {len(quality)}"
        f.write(f"@SEQ{i}\n")
        f.write(f"{seq}\n")
        f.write("+\n")
        f.write(f"{quality}\n")

print("FASTQ files created successfully")
print(f"Sequence length: {len(sequences_r1[0])}")
print(f"Quality length: {len(quality)}")
