include { SPLICEQ } from './modules/spliceq.nf'

params.projectDir = "/home/xiang/Projects/RNAseqAna"
params.gtf = "$projectDir/gencode.v38.primary_assembly.annotation.gtf"

workflow SPLICE_Q {
    // Read BAM files into a channel
    Channel
        .fromPath('aligned/*.bam')
        .map { file -> 
        def name = file.getName().replace('.bam', '')
        tuple(file, name)
        }
        .set { bam_ch }

    // Use the same GTF file for all BAM files
    SPLICEQ(bam_ch, params.gtf)
}

workflow {
    SPLICE_Q()
}