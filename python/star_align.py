import os
import subprocess
import re


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

alignment_folder = "aligned"
list = get_para(configfile = 'config_ORF52_K9_CLIP.txt')

fastq_files = get_fastq(configfile = 'config_ORF52_K9_CLIP.txt')

sample_number = int(len(fastq_files)/2)

outdir = list[3]
align_path = os.path.join(outdir, alignment_folder)
print(align_path)
#os.mkdir(align_path)
print(list)

for i in range(sample_number):
    subprocess.run(["./bash/align.sh", list[0], list[1], list[2], fastq_files[i*2], fastq_files[i*2+1], list[3], align_path])
