import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor

def process_chromosome(samtools_dir, bamfile_path, output_dir, chrom):
    output_bam = os.path.join(output_dir, f"{chrom}.bam")
    output_txt = os.path.join(output_dir, f"{chrom}.txt")
    try:
        subprocess.run([samtools_dir, 'view', '-b', bamfile_path, chrom, '-o', output_bam], check=True)
        with open(output_txt, 'w') as txt_file:
            subprocess.run([samtools_dir, 'pileup', output_bam], stdout=txt_file, check=True)
        # Use shell=True for the command that requires piping
        cmd = f"cut -f 2,4 {output_txt} > {output_txt}"
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    except subprocess.CalledProcessError:
        print(f"Read coverage file could not be generated for {chrom}")

def SamtoText(input_dir, current, bamfile_name, chromosomes):
    input_dir = os.path.abspath(input_dir)
    current = os.path.abspath(current)
    bamfile_path = os.path.join(current, input_dir, bamfile_name)
    
    output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
    last_dir = os.path.basename(os.path.normpath(bamfile_name)).split('.bam')[0]
    os.makedirs(last_dir, exist_ok=True)

    samtools_dir = "/usr/bin/samtools-0.1.8/samtools"

    try:
        subprocess.run([samtools_dir, 'index', bamfile_path], check=True)
    except subprocess.CalledProcessError:
        print("Index file could not be generated")
        return

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_chromosome, samtools_dir, bamfile_path, output_dir, chrom) for chrom in chromosomes]
        for future in futures:
            future.result()  # This will re-raise any exception that occurred during the execution

    print("Read coverage files generated for", bamfile_name)
