process BAM_PREPROCESSING {
    tag "$bam.simpleName"
    label 'process_medium'
    container 'biocontainers/picard:2.27.5'
    publishDir "${params.outdir}/bam_preprocessing/${bam.simpleName}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("_final.bam")) "bam/$filename"
            else if (filename.endsWith(".txt")) "metrics/$filename"
            else if (filename.endsWith(".pdf")) "metrics/$filename"
            else null
        }

    input:
    path bam
    path genome_fasta

    output:
    tuple val(bam.simpleName), path("*_final.bam"), emit: processed_bam
    path "*_metrics.txt", emit: metrics
    path "*_quality_metrics.txt", emit: quality_metrics
    path "*_insert_metrics.txt", emit: insert_metrics
    path "*_alignment_metrics.txt", emit: alignment_metrics
    path "*_gc_metrics.txt", emit: gc_metrics
    path "*.pdf", emit: plots

    script:
    """
    # Sort BAM file using Picard
    java -Xmx4g -jar /usr/picard/picard.jar SortSam \
        -I ${bam} \
        -O ${bam.simpleName}_sorted.bam \
        -SO coordinate \
        -CREATE_INDEX true \
        -VALIDATION_STRINGENCY LENIENT \
        -TMP_DIR .

    # Filter BAM for mapping quality using Picard
    java -Xmx4g -jar /usr/picard/picard.jar FilterSamReads \
        -I ${bam.simpleName}_sorted.bam \
        -O ${bam.simpleName}_filtered.bam \
        --FILTER includeAligned \
        --min-mapping-quality 20 \
        -CREATE_INDEX true \
        -VALIDATION_STRINGENCY LENIENT \
        -TMP_DIR .

    # Mark and remove duplicates using Picard
    java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates \
        -I ${bam.simpleName}_filtered.bam \
        -O ${bam.simpleName}_final.bam \
        -M ${bam.simpleName}_metrics.txt \
        -REMOVE_DUPLICATES true \
        -CREATE_INDEX true \
        -VALIDATION_STRINGENCY LENIENT \
        -TMP_DIR .

    # Collect quality metrics
    java -Xmx4g -jar /usr/picard/picard.jar QualityScoreDistribution \
        -I ${bam.simpleName}_final.bam \
        -O ${bam.simpleName}_quality_metrics.txt \
        -CHART ${bam.simpleName}_quality_metrics.pdf \
        -VALIDATION_STRINGENCY LENIENT

    # Collect insert size metrics (if paired-end)
    java -Xmx4g -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
        -I ${bam.simpleName}_final.bam \
        -O ${bam.simpleName}_insert_metrics.txt \
        -H ${bam.simpleName}_insert_metrics.pdf \
        -VALIDATION_STRINGENCY LENIENT \
        -ASSUME_SORTED true

    # Collect alignment metrics
    java -Xmx4g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
        -R ${genome_fasta} \
        -I ${bam.simpleName}_final.bam \
        -O ${bam.simpleName}_alignment_metrics.txt \
        -VALIDATION_STRINGENCY LENIENT

    # Collect GC bias metrics
    java -Xmx4g -jar /usr/picard/picard.jar CollectGcBiasMetrics \
        -I ${bam.simpleName}_final.bam \
        -O ${bam.simpleName}_gc_metrics.txt \
        -CHART ${bam.simpleName}_gc_metrics.pdf \
        -S ${bam.simpleName}_gc_summary.txt \
        -R ${genome_fasta} \
        -VALIDATION_STRINGENCY LENIENT

    # Clean up intermediate files
    rm ${bam.simpleName}_sorted.bam ${bam.simpleName}_filtered.bam
    """
} 