#!/usr/bin/env nextflow

include { STAR_INDEX } from './modules/star_align.nf'
include { ALIGN } from './modules/star_align.nf'
include { FASTQC } from './modules/fastqc.nf'
include { TRIM } from './modules/fastp_trim.nf'
include { FEATURECOUNT } from './modules/featurecount.nf'
include { MULTIQC } from './modules/multiqc.nf'

params.cpus = 12 
params.ram = 60000000000 /* ~60GB */

params.starindex = "$projectDir/star_index"
params.trimmeddir = "$projectDir/trimmed"

params.aligneddir = "$projectDir/aligned"
/* params.refgenome = "$projectDir/GRCh38.primary_assembly.genome.fa" */
params.gtf = "$projectDir/gencode.vM25.primary_assembly.annotation.gtf"

params.refgenomelink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz"
params.refgtflink = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz"

params.beddir = "$projectDir/annotation_bedfiles"

// Debugging statements to print parameter values
println "Project Directory: $projectDir"
println "CPUs: ${params.cpus}"
println "RAM: ${params.ram}"
println "STAR Index: ${params.starindex}"
println "Trimmed Directory: ${params.trimmeddir}"
println "Aligned Directory: ${params.aligneddir}"
println "GTF File: ${params.gtf}"


workflow RNASEQ {
	refgenome = file("${projectDir}/*.genome.fa")	

	
	Channel
	   	.fromFilePairs("${projectDir}/data/*{1,2}_001.fastq.gz", checkIfExists: true)
	   	.set { read_pairs_ch }
	read_pairs_ch.view()

	STAR_INDEX(refgenome)
	STAR_INDEX.out.view()


	FASTQC(read_pairs_ch)
	FASTQC.out.view()
	
	TRIM(read_pairs_ch)
	TRIM.out.view()

	ALIGN(STAR_INDEX.out.collect(), TRIM.out)
	ALIGN.out.view()

	FEATURECOUNT(params.gtf, ALIGN.out.collect())

	FEATURECOUNT.out.view()

	
    MULTIQC(".")
    MULTIQC.out.view()

}

workflow  {
	RNASEQ()
}
