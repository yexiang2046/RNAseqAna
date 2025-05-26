#!/usr/bin/env nextflow

process BEDTOOLS_GETFASTA {
    tag "Extract RNA sequences from ${bed_file.name}"
    
    input:
    path genome_fasta
    path bed_file
    val strand = true
    val output_name
    
    output:
    path "${output_name}.fasta", emit: fasta
    
    script:
    strand_opt = strand ? "-s" : ""
    """
    bedtools getfasta \\
        -fi ${genome_fasta} \\
        -bed ${bed_file} \\
        ${strand_opt} \\
        -split \\
        -name \\
        -fo ${output_name}.fasta
    """
} 