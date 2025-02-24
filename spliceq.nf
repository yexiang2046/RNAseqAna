include { SPLICEQ } from './modules/spliceq.nf'

workflow SPLICE_Q {

    // Read the GTF file into a single channel
    Channel.fromPath('*.gtf').set { annotation }

    // Read BAM files into a channel
    Channel.fromPath('*.bam').set { bamfile }

    // Use the same GTF file for all BAM files
    SPLICEQ(bamfile, annotation)
}

workflow {
    SPLICE_Q()
}