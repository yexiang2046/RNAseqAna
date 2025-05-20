#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BAM_PREPROCESSING } from '../modules/bam_preprocessing'
include { PIRANHA_PEAK_CALLING } from '../modules/piranha_peak_calling'


// Parameters
params.gtf = 'gencode.v38.primary_assembly.annotation.gtf'
params.bam_dir = null
params.rmsk = null
params.outdir = 'results'
params.piranha_params = ''
params.genome_fasta = null
params.viral_chr = 'NC_006998.1'

// Print usage message if required parameters are missing
if (!params.gtf || !params.bam_dir || !params.rmsk || !params.genome_fasta) {
    log.error """
    Required parameters missing!
    
    Usage:
    nextflow run main.nf --gtf GTF_FILE --bam_dir BAM_DIR --rmsk RMSK_BED --genome_fasta FASTA [options]
    
    Required arguments:
      --gtf          GTF file (e.g., gencode.v43.annotation.gtf)
      --bam_dir      Directory containing BAM files
      --rmsk         RepeatMasker BED file from UCSC
      --genome_fasta Reference genome FASTA file
    
    Optional arguments:
      --outdir            Output directory (default: 'results')
      --piranha_params    Additional parameters for Piranha
    """
    exit 1
}

// Process to extract genomic features from GTF
process EXTRACT_FEATURES {
    publishDir "${params.outdir}/features_bed", mode: 'copy'
    
    input:
    path gtf
    
    output:
    path "*.bed", emit: feature_beds
    
    script:
    """
    ${baseDir}/../bin/gtf_features_bedtools.sh -i $gtf -o .
    """
}



// Process to annotate peaks with genomic features
process ANNOTATE_FEATURES {
    publishDir "${params.outdir}/peak_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path feature_dir
    
    output:
    tuple val(sample_id), path("*_annotated_peaks.bed"), path("*_feature_overlap_summary.csv"), path("*_feature_overlap_plot.pdf")
    
    script:
    """
    # Make script executable and copy it to current directory
    cp ${projectDir}/../bin/annotate_peaks.sh .
    chmod +x annotate_peaks.sh
    
    # Run script from current directory
    ./annotate_peaks.sh \\
        -p $peaks \\
        -f $feature_dir \\
        -o .
    """
}

// Process to annotate peaks with repeat elements
process ANNOTATE_REPEATS {
    publishDir "${params.outdir}/repeat_annotations/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path rmsk
    
    output:
    tuple val(sample_id), path("*_rmsk_counts.bed"), path("*_rmsk_summary.csv"), path("*_rmsk_distribution.pdf")
    
    script:
    """
    # Make script executable and copy it to current directory
    cp ${projectDir}/bin/annotate_peaks_rmsk.sh .
    chmod +x annotate_peaks_rmsk.sh
    
    # Run script from current directory
    ./annotate_peaks_rmsk.sh \\
        -p $peaks \\
        -r $rmsk \\
        -o .
    """
}

// Process to count reads on peaks and calculate viral coverage
process COUNT_PEAK_READS {
    tag "$bam_meta"
    label 'process_high'
    publishDir "${params.outdir}/peak_counts/${bam_meta}", mode: 'copy'
    
    input:
    tuple val(bam_meta), path(bam)
    tuple val(bam_meta), path(peaks)
    val viral_chr
    
    output:
    tuple val(bam_meta), path("*_peak_counts.txt"), emit: peak_counts
    tuple val(bam_meta), path("*_viral_coverage.txt"), emit: viral_coverage
    tuple val(bam_meta), path("*_viral_depth.txt"), emit: viral_depth
    
    script:
    """
    # Count reads in peaks using bedtools
    bedtools coverage -a $peaks -b $bam -counts -s > ${bam_meta}_peak_counts.txt
    
    # Calculate viral genome coverage
    # First get viral genome length
    viral_length=\$(samtools view -H $bam | grep "SN:${viral_chr}" | awk '{print \$3}' | sed 's/LN://')
    
    # Calculate coverage on viral genome
    samtools index -@ 10 $bam
    samtools depth -r ${viral_chr} $bam > ${bam_meta}_viral_depth.txt
    
    # Calculate coverage statistics
    awk -v chr="${viral_chr}" -v len=\$viral_length '
    BEGIN {
        total_bases = 0
        covered_bases = 0
        total_depth = 0
        max_depth = 0
    }
    {
        if (\$3 > 0) {
            covered_bases++
            total_depth += \$3
            if (\$3 > max_depth) max_depth = \$3
        }
        total_bases++
    }
    END {
        avg_depth = (total_bases > 0) ? total_depth/total_bases : 0
        coverage = (total_bases > 0) ? (covered_bases/total_bases)*100 : 0
        print "Viral Genome Coverage Statistics"
        print "================================"
        print "Chromosome: " chr
        print "Genome length: " len " bp"
        print "Total bases sequenced: " total_bases
        print "Bases with coverage: " covered_bases
        print "Average depth: " avg_depth
        print "Maximum depth: " max_depth
        print "Coverage percentage: " coverage "%"
        print "Total reads: " total_depth
    }' ${bam_meta}_viral_depth.txt > ${bam_meta}_viral_coverage.txt
    """
}

// Main workflow
workflow {
    // Input channel for BAM files
    bam_ch = Channel.fromPath("${params.bam_dir}/*.bam")
    // genome_fasta = Channel.fromPath(params.genome_fasta)

    // Run the workflow
    BAM_PREPROCESSING(bam_ch)
    PIRANHA_PEAK_CALLING(BAM_PREPROCESSING.out.processed_bam)
    
    // Count reads on peaks and calculate viral coverage
    COUNT_PEAK_READS(
        BAM_PREPROCESSING.out.processed_bam,
        PIRANHA_PEAK_CALLING.out.peaks_bed,
        params.viral_chr
    )
    
    // Extract features from GTF
    // EXTRACT_FEATURES(gtf_file)
    
    // Annotate peaks with features
    // ANNOTATE_FEATURES(
    //     PIRANHA_PEAK_CALLING.out.peaks,
    //     EXTRACT_FEATURES.out.feature_beds.collect()
    // )
    
    // Annotate peaks with repeats
    // ANNOTATE_REPEATS(
    //     PIRANHA_PEAK_CALLING.out.peaks,
    //     rmsk_file
    // )
    
    // Generate final report
    // GENERATE_REPORT(
    //     ANNOTATE_FEATURES.out.collect()
    //         .mix(ANNOTATE_REPEATS.out.collect())
    //         .collect()
    // )
}