import re
import os
import subprocess
from glob import glob

def get_para(configfile):

    with open(configfile) as f:

        for line in f:
            line.rstrip()
            m = re.match(r"(Threads):(\d+)", line)
            idx = re.match(r"^(Index):(.*)$", line)
            dir = re.match(r"^(Working_dir):(.*)$", line)
            outdir = re.match(r"(Output_dir):(.*)$", line)
            if m:
                thread_num = m.group(2)
            elif idx:
                align_idx = idx.group(2)
            elif dir:
                work_dir = dir.group(2)
            elif outdir:
                out_dir = outdir.group(2)

    return thread_num, align_idx, work_dir, out_dir


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

# filter out mito DNA alignments
def filter_bam(bamfiles, thread, outdir):

    for f in bamfiles:
        sample = subprocess.check_output(["basename", f, ".out.bam"])
        dsample = sample.decode("utf8")
        #print(dsample)
        output_file = dsample.rstrip()+".filtered.bam"
        subprocess.run(["samtools", "index", "-@", thread, f])
        output_path = os.path.join(outdir, output_file)
        subprocess.run(["samtools", "view", "-b", "-h", "-@", thread, "-o", output_path, f, "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrGQ994935"])




def count_features(bamfiles, thread, outdir, outfile):
    bam_files = bamfiles
    subprocess.run(["featureCounts", "-p", "--countReadPairs", "-t", "exon", "-g", "gene_id", "-a", "annotation.gtf", "-o", os.path.join(outdir, outfile), bam_files])



list = get_para(configfile = 'config_ORF52_K9_CLIP.txt')
Threads = list[0]
Work_dir = list[2]
Output_dir = list[3]
Alignment_dir = os.path.join(Output_dir, "aligned")
Alignment_files = glob(os.path.join(Alignment_dir, "*.out.bam"))

#filter_bam(bamfiles = Alignment_files, thread=Threads, outdir=Alignment_dir)
Filtered_files = glob(os.path.join(Alignment_dir, "*.filtered.bam"))
print(Filtered_files)
Outfile = "featureCounts.txt"
for f in Filtered_files:
    count_features(bamfiles = Filtered_files, thread=Threads, outdir=Output_dir, outfile=f)
print("Please find output at:")
print(Output_dir)
print("Alignment results:")
print(Alignment_files)
