
include { SPLICEQ } from './modules/spliceq.nf'
workflow SPLICE_Q {

    // Read the GTF file into a single channel
    Channel.fromPath('*.gtf').set { annotation }

    // Read BAM files into a channel
    Channel.fromPath('*.bam').set { bamfile }

    // Use the same GTF file for all BAM files
    bamfile
        .combine(annotation)
        .set { bamfile_with_annotation }

    SPLICEQ(bamfile_with_annotation)
}

workflow {
    SPLICE_Q()
}