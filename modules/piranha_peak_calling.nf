process PIRANHA_PEAK_CALLING {
    tag "$bam_meta"
    label 'process_high'
    publishDir "${params.outdir}/piranha_output/${bam_meta}", 
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith("_peaks.bed")) "peaks/$filename"
            else if (filename.endsWith("_peaks.tsv")) "peaks/$filename"
            else if (filename.endsWith("_converted.bed")) "bed/$filename"
            else if (filename.endsWith("_piranha_metrics.txt")) "metrics/$filename"
            else null
        }
    
    input:
    tuple val(bam_meta), path(processed_bam)
    
    output:
    tuple val(bam_meta), path("*_peaks.tsv"), emit: peaks
    tuple val(bam_meta), path("*_peaks.bed"), emit: peaks_bed
    path "*_piranha_metrics.txt", emit: metrics
    
    script:
    def piranha_opts = params.piranha_params ?: '-s -b 20 -d ZeroTruncatedNegativeBinomial -p 0.01 -u 100'
    """
    # Run Piranha on the processed BAM file
    Piranha $piranha_opts ${processed_bam} > ${bam_meta}_peaks.tsv
    
    # Convert Piranha output to BED format with proper scoring and p-value in name
    awk -v OFS='\t' '
        NR>1 {
            chr=\$1
            start=\$2
            end=\$3
            pval=\$7
            # Format name as sample_peak_number_pvalue
            name=sprintf("%s_peak_%d_pval%.2e", "'${bam_meta}'", NR-1, pval)
            score=\$4
            strand=\$6
            print chr, start, end, name, score, strand
        }' ${bam_meta}_peaks.tsv > ${bam_meta}_peaks.bed
    
    # Generate peak calling metrics
    echo "Peak Calling Metrics for ${bam_meta}" > ${bam_meta}_piranha_metrics.txt
    echo "--------------------------------" >> ${bam_meta}_piranha_metrics.txt
    echo "Total peaks called: \$(wc -l < ${bam_meta}_peaks.bed)" >> ${bam_meta}_piranha_metrics.txt
    echo "Parameters used: $piranha_opts" >> ${bam_meta}_piranha_metrics.txt
    echo "" >> ${bam_meta}_piranha_metrics.txt
    echo "Peak Score Distribution:" >> ${bam_meta}_piranha_metrics.txt
    awk '
        BEGIN {print "Min\tQ1\tMedian\tQ3\tMax"}
        {
            print \$1
        }
        END {
            # Use sort command to get sorted values
            cmd = "sort -n"
            while ((cmd | getline line) > 0) {
                scores[++n] = line
            }
            close(cmd)
            print scores[1] "\t" \
                  scores[int(n*0.25)] "\t" \
                  scores[int(n*0.5)] "\t" \
                  scores[int(n*0.75)] "\t" \
                  scores[n]
        }
    ' ${bam_meta}_peaks.bed >> ${bam_meta}_piranha_metrics.txt
    """
} 