include { SPLICEQ } from './modules/spliceq.nf'

// Define default parameters
params.bam_dir = 'aligned'
params.gtf = null
params.outdir = 'results'

workflow SPLICE_Q {
    // Parameter validation
    if (!params.gtf) {
        error "GTF file not specified! Please provide --gtf parameter."
    }

    // Read BAM files into a channel
    Channel
        .fromPath("${params.bam_dir}/*.bam")
        .ifEmpty { error "No BAM files found in ${params.bam_dir}" }
        .map { file -> 
            def name = file.getName().replace('.bam', '')
            tuple(name, file)
        }
        .set { bam_ch }

    // Log input parameters
    log.info """
    ==============================================
    SPLICE-Q Pipeline
    ==============================================
    BAM directory : ${params.bam_dir}
    GTF file     : ${params.gtf}
    Output dir   : ${params.outdir}
    """

    // Run SPLICE-Q
    SPLICEQ(bam_ch, file(params.gtf))
}

workflow {
    SPLICE_Q()
}