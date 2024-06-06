import os
import sys
import time
from multiprocessing import Pool

def process_chromosome(args):
    samtools_dir, input_dir, current, bamfile_name, chrom, output_dir = args
    try:
        cmd2 = f"{samtools_dir} view -b {os.path.join(current, input_dir, bamfile_name)} {chrom} -o {os.path.join(current, output_dir, chrom + '.bam')}"
        cmd3 = f"{samtools_dir} pileup {os.path.join(current, output_dir, chrom + '.bam')} | cut -f 2,4 > {os.path.join(current, output_dir, chrom + '.txt')}"
        command = f"{cmd2}; {cmd3}"
        os.system(command)
    except ValueError:
        print(f"Read coverage file for chromosome {chrom} could not be generated")
        sys.exit()

def SamtoText(input_dir, current, bamfile_name, chromosomes, cores):
    output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
    last_dir = os.path.basename(os.path.normpath(bamfile_name)).split('.bam')[0]
    os.makedirs(last_dir, exist_ok=True)

    samtools_dir = "/usr/bin/samtools-0.1.8/samtools"
    try:
        cmd1 = f"{samtools_dir} index {os.path.join(current, input_dir, bamfile_name)}"
        os.system(cmd1)
    except ValueError:
        print("Index file could not be generated")
        sys.exit()

    args_list = [(samtools_dir, input_dir, current, bamfile_name, chrom, output_dir) for chrom in chromosomes]

    with Pool(int(cores)) as pool:
        pool.map(process_chromosome, args_list)

    # Clean up intermediate files
    for chrom in chromosomes:
        bam_path = os.path.join(current, output_dir, chrom + '.bam')
        if os.path.exists(bam_path):
            os.remove(bam_path)

    return