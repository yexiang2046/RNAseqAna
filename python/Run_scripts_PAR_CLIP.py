# make sure to have cutadapt, samtools, and deeptools installed under the environment



import re
import subprocess
import os
from glob import glob
import itertools
#from bioinfodit.analysis import gff


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


    return thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_dir, fastq_dir, genome, annotation, featureCounts_out


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


def index_bam(bamfile, thread):

        subprocess.run(["samtools", "index", "-@", thread, bamfile])



def filter_bam(bamfile, thread, outdir, output_file, genome):
    if genome == "hg38":
        # output_path = os.path.join(outdir, output_file)
        # filtered on combined genome of human and HHV8
        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrGQ994935"])
    elif genome == "mm10":

        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_file, bamfile, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"])


def bam_to_bw(bamfile, thread, output_file):

    subprocess.run(["bamCoverage", "-b", bamfile, "-o", output_file, "--binSize", "10", "-p", thread, "--normalizeUsing", "CPM"])


def count_features(bamfiles, gtf, thread, count_outfile):
    bam_files = bamfiles
    subprocess.run(["featureCounts", "-T", thread, "-p", "-t", "exon", "-g", "gene_id", "-a", gtf, "-o", count_outfile]  + bam_files)





Para = get_para(configfile="config_PAR_CLIP.txt")
thread_num = Para[0]
align_idx = Para[1]
out_dir = Para[2]
trimmer_dir = Para[3]
alignment_dir = Para[4]
pair_end = Para[5]
trimmed_output = Para[6]
fastq_dir = Para[7]
genome = Para[8]
annotation = Para[9]
counts_outfile = Para[10]


######################################
# change the config file accorddingly
fastq_files = get_fastq(configfile="config_PAR_CLIP.txt")

samples = []

for f in fastq_files:
    samples.append(os.path.basename(f).split('_', 1)[0])

alignment_path = os.path.join(out_dir, alignment_dir)
trimmed_path = os.path.join(out_dir, trimmed_output)

if not os.path.exists(trimmed_path):
    os.mkdir(trimmed_path)

if not os.path.exists(alignment_path):
    os.mkdir(alignment_path)

bw_dir = os.path.join(out_dir, "bw_files")
if not os.path.exists(bw_dir):
    os.mkdir(bw_dir)


print(Para)
print(fastq_files)
if (pair_end == "Yes"):
    sample_num = int(len(fastq_files)/2)
    for i in range(sample_num):
    #    subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_path, fastq_files[i*2], fastq_files[i*2+1]])
        print(fastq_files[i*2], fastq_files[i*2+1])
else:
    for i in range(len(fastq_files)):
    #    subprocess.run(["./bash/step1_reads_trimming.sh", thread_num, align_idx, out_dir, trimmer_dir, alignment_dir, pair_end, trimmed_path, os.path.join(fastq_dir, fastq_files[i])])
        print(fastq_files[i])
    print(os.path.join(fastq_dir, fastq_files[1]))

sample_num = int(len(fastq_files)/2)
for i in range(sample_num):
    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
#    print(file)
    trimmed_file1 = "".join([file, "_S1_L005_R1_001_val_1.fq.gz"])
    trimmed_files1 = os.path.join(trimmed_path, trimmed_file1)
    trimmed_file2 = "".join([file, "_S1_L005_R2_001_val_2.fq.gz"])
    trimmed_files2 = os.path.join(trimmed_path, trimmed_file2)
#    print(trimmed_files1, trimmed_files2)
    #subprocess.run(["./bash/align.sh", thread_num, align_idx, out_dir, trimmed_files1, trimmed_files2, alignment_path])


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    bamfile = "".join([alignment_path, "/", file, "_S1_L005_R1_001_val_1.fq.gzAligned.sortedByCoord.out.bam"])
    outputfile = "".join([alignment_path, "/", file, ".filtered.bam"])
    #index_bam(bamfile = bamfile, thread = thread_num)
    #filter_bam(bamfile = bamfile, thread = thread_num, outdir = alignment_path, output_file = outputfile, genome = genome)

for i in range(sample_num):

        file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
        filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
        #index_bam(filtered_bam, thread = thread_num)


for i in range(sample_num):

    file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
    filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
    bigwig_file = "".join([bw_dir, "/", file, ".bw"])

    bam_to_bw(bamfile = filtered_bam, thread=thread_num, output_file=bigwig_file)

bam_list =[]
for i in range(sample_num):

        file = os.path.basename(fastq_files[i*2]).split('_', 1)[0]
        filtered_bam = "".join([alignment_path, "/", file, ".filtered.bam"])
        bam_list.append(filtered_bam)


GTF = annotation
GTF_file = os.path.join(out_dir, GTF)
count_file = counts_outfile
count_path_file = os.path.join(out_dir, count_file)
#print(count_file)
#print(count_path_file)
#print(bam_list)
count_features(bamfiles = bam_list, gtf=GTF_file, thread=thread_num, count_outfile=count_path_file)


## Generate report files
subprocess.run(["bash", "source", "/home/xiang/miniconda3/bin/activate"])
subprocess.run(["multiqc", out_dir, "-o", out_dir])
