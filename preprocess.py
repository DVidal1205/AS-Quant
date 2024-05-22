import os
import sys
import subprocess
import time

def SamtoText(input_dir, current, bamfile_name, chromosomes):
	output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
	last_dir = os.path.basename(os.path.normpath(bamfile_name)).split('.bam')[0]
	os.makedirs(last_dir, exist_ok=True)

	samtools_dir = "/usr/bin/samtools-0.1.8/samtools"
	try:
		cmd1 = [samtools_dir, "index", os.path.join(current, input_dir, bamfile_name)]  # make samtools index filename.bam.bai
		subprocess.run(cmd1, check=True)
		# print("bam index file generated.")
	except subprocess.CalledProcessError:
		print("Index file could not be generated")
		sys.exit()
		for chrom in chromosomes:
			# print("Chrom ", chrom, "... ...")
			tt = time.time()
			cmd2 = [samtools_dir, "view", "-b", os.path.join(current, input_dir, bamfile_name), chrom, "-o",
					os.path.join(current, output_dir, chrom + ".bam")]
			cmd3 = [samtools_dir, "pileup", os.path.join(current, output_dir, chrom + ".bam")]
			cmd4 = ["cut", "-f", "2,4"]
			cmd5 = [">", os.path.join(current, output_dir, chrom + ".txt")]
			command = cmd2 + [";", cmd3] + ["|"] + cmd4 + [">"] + cmd5
			# print(command)
			try:
				command_str = ' '.join([str(elem) for elem in command])
				command_str = command_str.replace('[', '')			
				command_str = command_str.replace(']', '')
				command_str = command_str.replace('\'', '')
				command_str = command_str.replace(',', '')

				print(command_str)
				subprocess.run(command_str, shell=True, check=True, executable="/bin/bash")
				print("Read coverage files generated for", bamfile_name)
			except subprocess.CalledProcessError:
				print("Read coverage file could not be generated")
				sys.exit()

	return
