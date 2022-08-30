import re
import subprocess
import os


def get_para(configfile):

    with open(configfile) as f:

        for line in f:
            line.rstrip()
            m = re.match(r"(Threads):(\d+)", line)
            idx = re.match(r"^(Index):(.*)$", line)
            outdir = re.match(r"(Output_dir):(.*)$", line)
            trim = re.match(r"(TrimGalore):(.*)$", line)
            align = re.match(r"(Alignment_dir):(.*)$", line)
            pairend = re.match(r"(PAIREND):(.*)$", line)
            trimmed = re.match(r"(Trimmed_dir):(.*)$", line)
            fastqdir = re.match(r"(Fastq_dir):(.*)", line)
            species = re.match(r"(GENOME):(.*)", line)
            gtf = re.match(r"(ANNOTATION):(.*)", line)
            counts = re.match(r"(FEATURECOUNTS_OUT):(.*)", line)
            aligner = re.match(r"(ALIGNER):(.*)", line)
            if m:
                thread_num = m.group(2)
            elif idx:
                align_idx = idx.group(2)
            elif outdir:
                out_dir = outdir.group(2)
            elif trim:
                trimmer_dir = trim.group(2)
            elif align:
                alignment_dir = align.group(2)
            elif pairend:
                pair_end = pairend.group(2)
            elif trimmed:
                trimmed_dir = trimmed.group(2)
            elif fastqdir:
                fastq_dir = fastqdir.group(2)
            elif species:
                genome = species.group(2)
            elif gtf:
                annotation = gtf.group(2)
            elif counts:
                featureCounts_out = counts.group(2)
            elif aligner:
                aligners = aligner.group(2)
    return thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_dir, fastq_dir, genome, annotation, featureCounts_out, aligners


def get_fastq(configfile):

    with open(configfile) as f:
        files = []
        for line in f:
            line.rstrip()
            fq = re.match(r"^(.*fastq\.gz$)", line)
            if fq:
                file = fq.group(1)
                files.append(file)

    return files

def generate_hisat2_idx(fasta_file, idx_folder):

    subprocess.run(["hisat2-build", fasta_file, idx_folder])

def hisat2_align(thread_num, idx, reads_1, reads_2, output_sam):

    subprocess.run(["hisat2", "--no-unal", "--rna-strandness", "FR", "-p", thread_num, "-x", idx, "-1", reads_1, "-2", reads_2, "-S", output_sam])


def sam2bam_q10(sam, bam, thread_num=8):
    subprocess.run(["samtools", "view", "-@", "{thread_num}", "-b", "-q", "10", "-o", bam, sam])

def sortbam(bam, sortedBam, thread_num):
    subprocess.run(["samtools", "sort", "-@", thread_num, "-o", sortedBam, bam])


def index_bam(bamfile, thread):

        subprocess.run(["samtools", "index", "-@", thread, bamfile])



def filter_bam(bamfile, thread, outdir, output_file, genome):
    if genome == "hg38_KSHV":
        # output_path = os.path.join(outdir, output_file)
        # filtered on combined genome of human and HHV8
        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrGQ994935"])
    elif genome == "mm10":

        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"])
    elif genome == "hg38":

        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"])
    
    elif genome == "hg38_VACV":

        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "NC_006998.1"])

def bam_to_bw(bamfile, thread, output_file, norm_method = "CPM"):

    subprocess.run(["bamCoverage", "-b", bamfile, "-o", output_file, "--binSize", "20", "-p", thread, "--normalizeUsing", norm_method])


def count_features(bamfiles, gtf, thread, count_outfile):
    bam_files = bamfiles
    subprocess.run(["featureCounts", "-T", thread, "-p", "-t", "exon", "-g", "gene_id", "-a", gtf, "-o", count_outfile]  + bam_files)
