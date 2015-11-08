
#!/usr/bin/env python

import commands
import ntpath
import sys
import os
import fileinput


help_menue = "\nHierarchical Genome Assembly version 1.0.0\n" + \
"--Prerequisite-- \n" + \
"Velvet should be installed with option \'LONGSEQUENCES=1\' in the make command, to allow Velvet to accept contigs (long sequences) as input; as well and \'MAXKMERLENGTH=111\' as the default max kmer length is 31\n\n" + \
"--Parameters--: \n" + \
" -velvet  path      Path to the directory where velvet bineries (velveth, velvetg) are located \n" + \
" -spades  path      Path to the directory where spades.py file is located\n" + \
" -PA      string    select one of {SPAdes, velvet} that will assemble the partitions\n" + \
" -P12     file      interlaced fastq file of paired reads that will be used in the partiton step, usually raw or clean (corrected) reads\n" + \
" -R12     file      interlaced fastq file of paired reads that will be used in the re-assembly step, usually raw or clean (corrected) reads\n" + \
" -ins     int       insert size of the fragments\n" + \
" -std     int       standard deviation of the inset size\n" + \
" -P       int       number of partitions\n" + \
" -Pkmer   int(odd)  kmer value for assembling the partitions\n" + \
" -Rkmer   int(odd)  kmer value for re-assembly step \n" + \
" -t       int       number of threads to be used for re-assembly step using SPAdes assembler default 1.\n" + \
" -out     Path      output path\n" + \
" -h                 print option menue\n" 


def splitting_fastq(reads_file, part_num):
		print 'splitting fastq file into ' + str(part_num) + " partitions"
		num_lines = 0
		for i in fileinput.input(reads_file):
			num_lines += 1

		part_size = num_lines/part_num
		part_size = part_size/8 * 8
		
		
		for i in range (1, part_num + 1):
			out1 = open( ntpath.basename(reads_file)[:-5] + "_part_" + str(i) +".fastq" , 'w')
			num_lines = 0
			for j in fileinput.input(reads_file):
				num_lines += 1
				if num_lines > (i - 1) * part_size and num_lines <= i * part_size:				
					out1.writelines(j)

def sanity_check(P_fastq_file, R_fastq_file, velvet, spades):
	fastq_flag = 0
	i = 0
	t = 0
	for line in fileinput.input(P_fastq_file):
		tmp = line.strip().split('\t')
		if i == 0 and tmp[0][0] == "@":
			fastq_flag += 1
		if i == 2 and tmp[0][0] == "+":
			fastq_flag += 1
		if i == 1:
			t = len(tmp[0])
		if i == 3 and len(tmp[0]) == t:
			fastq_flag += 1 
		if i == 4:
			break
		i += 1
	fileinput.close()
	if fastq_flag !=3 :
		print "\n Wrong fastq format for the partition step, please provide a fastq file with the standard format\n"
		exit()

	fastq_flag = 0
	i = 0
	t = 0
	for line in fileinput.input(R_fastq_file):
		tmp = line.strip().split('\t')
		if i == 0 and tmp[0][0] == "@":
			fastq_flag += 1
		if i == 2 and tmp[0][0] == "+":
			fastq_flag += 1
		if i == 1:
			t = len(tmp[0])
		if i == 3 and len(tmp[0]) == t:
			fastq_flag += 1 
		if i == 4:
			break
		i += 1
	fileinput.close()
	if fastq_flag !=3 :
		print "\n Wrong fastq format for the re-assembly step, please provide a fastq file with the standard format\n"
		exit()



	if "not found" in commands.getstatusoutput(velvet)[1]:
		print "\n Wrong velvet package; velveth can't be found in the given path\n"
		exit()

	if "not found" in commands.getstatusoutput(spades)[1]:
		print "\n Wrong SPAdes package; spade.py can't be found in the given path\n"
		exit()
	 


def run_command(parm):
	
	if ("-h" in parm or "-help" in parm or len(parm) == 0) :
		print help_menue
	
	else:
		velvet_path = parm[parm.index("-velvet") + 1]
		spades_path = parm[parm.index("-spades") + 1]
		partitions_fastq_file  = parm[parm.index("-P12") + 1]
		R_assembly_fastq_file  = parm[parm.index("-R12") + 1]


		selected_part_assem = parm[parm.index("-PA") + 1]
		if not (selected_part_assem == 'SPAdes' or selected_part_assem =='velvet'):
			print "\nPlease choose between SPAdes or velvet to assemble the partitions\n"
			exit()

		try:
			insert_size  = int(parm[parm.index("-ins") + 1])
		except ValueError:
			print "\nPlease provide integer insert size\n"
			exit()
		try:
			std  = int(parm[parm.index("-std") + 1])
		except ValueError:
			print "\nPlease provide integer standard deviation\n"
			exit()

		try:
			num_parts = int(parm[parm.index("-P") + 1])
		except ValueError:
			print "\nPlease provide integer number of partitions\n"
			exit()

		try:
			p_kmer = int(parm[parm.index("-Pkmer") + 1])
		except ValueError:
			print "\nPlease provide integer kmer size for partitions assembly\n"
			exit()

		try:
			r_kmer = int((parm[parm.index("-Rkmer") + 1]))
		except ValueError:
			print "\nPlease provide integer kmer size for re-assembly process\n"
			exit()
		try:
			tt = parm.index("-t")
			try:
 				t = int((parm[parm.index("-t") + 1]))
				t = parm[parm.index("-t") + 1]
			except ValueError:
				print "\nPlease provide integer number for SPAdes threading\n"
				exit()

		except ValueError:
			t = "1"

		insert_size  = parm[parm.index("-ins") + 1]
		std  = parm[parm.index("-std") + 1]
		num_parts = int(parm[parm.index("-P") + 1])
		p_kmer = parm[parm.index("-Pkmer") + 1]
		r_kmer = parm[parm.index("-Rkmer") + 1]
		out= parm[parm.index("-out") + 1]
		

# Check the sanity of fastq file, spade command, and velvet command.
		sanity_check(partitions_fastq_file, R_assembly_fastq_file, velvet_path+"/velveth", spades_path+"/spades.py")

		os.system("mkdir " + out)
		os.chdir(out)
		out_path = os.popen("pwd").read().split()[0]
		os.system("cd " + out_path)

# Split the fastq file into the given partitions
		splitting_fastq(partitions_fastq_file, num_parts)

		

		if selected_part_assem == "SPAdes":
			for p in range(1, num_parts + 1):
				os.chdir( out_path)
				os.system("mkdir part_" + str(p) + "_assembly")
				os.chdir("part_" + str(p) + "_assembly")
				
				print "Partition " + str(p) + " assembly started"
				os.popen(spades_path + "/spades.py -t " + t + " -k " + p_kmer + " --12 " + out_path + "/" + ntpath.basename(partitions_fastq_file)[:-5] + "_part_" + str(p) +".fastq" + " -o ./ ").read()

		
		if selected_part_assem == "velvet":
			for p in range(1, num_parts + 1):
				os.chdir( out_path)
				os.system("mkdir part_" + str(p) + "_assembly")
				os.chdir("part_" + str(p) + "_assembly")
				
				print "Partition " + str(p) + " assembly started"

				os.popen(velvet_path+"/velveth ./ " + p_kmer + " -fastq -shortPaired "+ out_path + "/" + ntpath.basename(partitions_fastq_file)[:-5] + "_part_" + str(p) +".fastq").read()
				os.popen(velvet_path+"/velvetg ./ -exp_cov auto -ins_length "+ insert_size +" -ins_length_sd " + std + " -scaffolding no").read()
				os.system("mv contigs.fa contigs.fasta")







		print "Merging contigs from all parts"
		os.chdir( out_path)
		os.system("mkdir merged_contigs")
		os.chdir("merged_contigs")

		out = open( "merged_contigs.fasta" , 'w')
		for p in range(1, num_parts + 1):
			for line in fileinput.input(out_path+"/"+"part_" + str(p) + "_assembly/contigs.fasta"):
				out.writelines(line)
		out.close

		print "Combining contigs using velvet, 31 as kmer"

		os.chdir( out_path)
		os.system("mkdir combined_contigs")
		os.chdir("combined_contigs")
		
		os.popen(velvet_path+"/velveth ./ 31 -fasta -long "+ out_path +"/merged_contigs/merged_contigs.fasta").read()
		os.popen(velvet_path+"/velvetg ./ -exp_cov "+ str(num_parts) + " -scaffolding no").read()
		os.system("mv contigs.fa combined_contigs.fasta")

		print "Re-assembling reads with merged contigs(HGA(merged))"
		os.chdir( out_path)
		os.system("mkdir HGA_merged ")
		os.chdir("HGA_merged")
		os.popen(spades_path + "/spades.py -t " + t + " -k " + r_kmer + " --12 " + R_assembly_fastq_file  + " --trusted-contigs " + out_path+"/merged_contigs/merged_contigs.fasta  -o ./ ").read()

		print "Re-assembling reads with combined contigs(HGA(combined))"
		os.chdir( out_path)
		os.system("mkdir HGA_combined")
		os.chdir("HGA_combined")
		os.popen(spades_path + "/spades.py -t " + t + " -k " + r_kmer + " --12 " + R_assembly_fastq_file  + " --trusted-contigs " + out_path+"/combined_contigs/combined_contigs.fasta  -o ./ ").read()


	
		
run_command(sys.argv[1:])







