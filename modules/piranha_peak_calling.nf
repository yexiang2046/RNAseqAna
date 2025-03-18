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
    path "*_unique.bed", emit: unique_bed
    
    script:
    def piranha_opts = params.piranha_params ?: '-s -b 20 -d ZeroTruncatedNegativeBinomial -p 0.00001 -u 100'
    """
    # Convert BAM to BED
    bedtools bamtobed -i $bam > ${bam.simpleName}_converted.bed
    
    # Filter out reads not uniquely mapped (using single quotes for awk)
    awk '\$5 == "255" {print \$0}' ${bam.simpleName}_converted.bed > ${bam.simpleName}_unique.bed
    
    # Run Piranha on the sorted BED file
    Piranha $piranha_opts ${bam.simpleName}_unique.bed > ${bam.simpleName}_peaks.tsv
    
    # Clean up intermediate files
    rm ${bam.simpleName}_converted.bed
    """
} 