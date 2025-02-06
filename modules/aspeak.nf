process ASPEAK {
    debug true
    container 'biocontainers/aspeak:v1.0.0'
    publishDir "${projectDir}/aspeak_out", mode:'copy'

    # take rnaseq_bamfile, rip_bamfile, and beddir as input
    input:
    tuple path(rnaseq_bamfile), path(rip_bamfile), path(beddir)

    # output to aspeak_out directory
    output:
    path "*.txt"


    script:
    """
    ./ASPeak/scripts/aspeak.pl -lib ${rnaseq_bamfile} -lib ${rip_bamfile} -beddir ${beddir} -outdir ${projectDir}/aspeak_out -rnaseq ${rnaseq_bamfile}
    """
}