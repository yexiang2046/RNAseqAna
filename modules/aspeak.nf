process ASPEAK {
    debug true
    publishDir "${projectDir}/aspeak_out", mode:'copy'

    input:
    tuple path(rnaseq_bamfile), path(rip_bamfile), path(beddir)

    output:
    path "*.txt"


    script:
    """
    ./ASPeak/scripts/aspeak.pl -lib ${rnaseq_bamfile} -lib ${rip_bamfile} -beddir ${beddir} -outdir ${projectDir}/aspeak_out -rnaseq ${rnaseq_bamfile}
    """
}