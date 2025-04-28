#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples = "samples.txt"
// sample sheet
// sample_id, treatment_bam, control_bam

params.outdir  = "macs2_results"
params.genome  = "hs" // or "mm", or provide effective genome size

process MACS2_PEAKCALL {
    tag "${sample_id}"
    container "quay.io/biocontainers/macs2:2.2.9.1--py311haab0aaa_3"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(treat_bam), path(ctrl_bam)

    output:
    path("${sample_id}_macs2_peaks.narrowPeak")
    path("${sample_id}_macs2_peaks.xls")
    path("${sample_id}_macs2_summits.bed")

    script:
    """
    macs2 callpeak \
        -t ${treat_bam} \
        -c ${ctrl_bam} \
        -f BAM \
        -g ${params.genome} \
        -n ${sample_id}_macs2 \
        --outdir . \
        --keep-dup all \
        -q 0.05
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.samples)
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
            tuple(row.sample, file(row.treatment_bam), file(row.control_bam))
        }
        | MACS2_PEAKCALL
}