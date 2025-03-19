process BAM_PREPROCESSING {
    tag "$bam.simpleName"
    label 'process_medium'
    container 'mgibio/samtools:1.9'
    publishDir "${params.outdir}/bam_preprocessing/${bam.simpleName}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("_final.bam")) "bam/$filename"
            else if (filename.endsWith(".txt")) "metrics/$filename"
            else null
        }

    input:
    path bam

    output:
    tuple val(bam.simpleName), path("*_final.bam"), emit: processed_bam
    path "*_stats.txt", emit: stats
    path "*_flagstat.txt", emit: flagstat
    path "*_idxstats.txt", emit: idxstats

    script:
    """
    # Sort BAM file
    samtools sort -@ 8 -m 4G \
        -o ${bam.simpleName}_sorted.bam \
        ${bam}

    # Index sorted BAM
    samtools index ${bam.simpleName}_sorted.bam

    # Filter for mapping quality and proper pairs
    samtools view -@ 8 -b -q 255 -F 0x4 \
        ${bam.simpleName}_sorted.bam \
        > ${bam.simpleName}_filtered.bam

    # Mark and remove duplicates
    samtools markdup -@ 8 -r \
        ${bam.simpleName}_filtered.bam \
        ${bam.simpleName}_final.bam

    # Index final BAM
    samtools index ${bam.simpleName}_final.bam

    # Collect various statistics
    samtools stats ${bam.simpleName}_final.bam > ${bam.simpleName}_stats.txt
    samtools flagstat ${bam.simpleName}_final.bam > ${bam.simpleName}_flagstat.txt
    samtools idxstats ${bam.simpleName}_final.bam > ${bam.simpleName}_idxstats.txt

    # Clean up intermediate files
    rm ${bam.simpleName}_sorted.bam ${bam.simpleName}_sorted.bam.bai ${bam.simpleName}_filtered.bam
    """
} 