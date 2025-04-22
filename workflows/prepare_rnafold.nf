include { BEDTOOLS_GETFASTA } from './modules/bedtools_getfasta'

workflow {
    BEDTOOLS_GETFASTA(
        genome_fasta,
        bed_file,
        true,  // strand (optional, defaults to true)
        "extracted_rna_sequences"  // output_name
    )
}