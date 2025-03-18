process PIRANHA_PEAK_CALLING {
    tag "$bam.simpleName"
    label 'process_high'
    container 'xiang2019/piranha:v1.0.0'
    publishDir "${params.outdir}/piranha_output/${bam.simpleName}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("_peaks.bed")) "peaks/$filename"
            else if (filename.endsWith("_peaks.tsv")) "peaks/$filename"
            else if (filename.endsWith("_converted.bed")) "bed/$filename"
            else null
        }
    
    input:
    path bam
    
    output:
    tuple val(bam.simpleName), path("*_peaks.tsv"), emit: peaks
    
    script:
    def piranha_opts = params.piranha_params ?: '-s -b 20 -d ZeroTruncatedNegativeBinomial -p 0.01 -u 100'
    """
    # Sort bam file by coordinates
    samtools sort -@ 8 -o ${bam.simpleName}_sorted.bam ${bam}
    
    # Filter for uniquely mapped reads (MAPQ >= 20)
    samtools view -@ 8 -b -q 20 ${bam.simpleName}_sorted.bam > ${bam.simpleName}_unique.bam
    
    # Mark and remove duplicates
    samtools markdup -@ 8 -r ${bam.simpleName}_unique.bam ${bam.simpleName}_filtered.bam
    
    # Index the filtered BAM
    samtools index ${bam.simpleName}_filtered.bam
    
    # Run Piranha on the filtered BAM file
    Piranha $piranha_opts ${bam.simpleName}_filtered.bam > ${bam.simpleName}_peaks.tsv
    
    # Clean up intermediate files
    rm ${bam.simpleName}_sorted.bam ${bam.simpleName}_unique.bam ${bam.simpleName}_filtered.bam ${bam.simpleName}_filtered.bam.bai
    """
} 