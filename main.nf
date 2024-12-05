#!/usr/bin/env nextflow

params.cpus = 8



params.fqdir = "$projectDir/data"
params.trimeddir = "$projectDir/trimmed"
params.aligneddir = "$projectDir/aligned"
params.refgenome = "$projectDIr/GRCh38.primary.genome.fa"

process trim{
	input:
	path	fqdir

	output:
	path	trimmeddir

	script:
	"""
	./bash/fastp_trim.sh fqdir _S1_L005_R1_001.fastq.gz _S1_L005_R2_001.fastq.gz trimmeddir
	"""
}



process align{
	input:
	path	trimmeddir

	output:
	path	aligneddir

	script:
	"""
	./bash/star_align.sh trimmeddir .R1.fastp.fastq.gz .R2.fastp.fastq.gz STAR_IDX aligneddir refgenome
	"""
}


workflow {
	trim()
	align()
}
