#script to use extract_features_from_gtf.py to extract features from gtf file
#!/bin/bash

# Extract features from gtf file
extract-transcript-regions/extract_transcript_regions.py \
    -i gencode.v38.primary_assembly.annotation.gtf \
    -o gencode.v38 \
    --gtf