# To change this template, choose Tools | Templates
# and open the template in the editor.
import random
import string
import heapq
import math
from collections import defaultdict
from datetime import datetime
import sys
import os


reads_array = []
compressed_genome = ""
uncompressed_genome = ""
kmer_dict = defaultdict(int)
kmer_sorted_freq_dict = defaultdict(list)
bases_freq_dict = defaultdict(int)
None_DNA_bases_freq_dict = defaultdict(int)
RLE_numbers_dict = defaultdict(int)


code_bits = []
bits_counter = 0
genome_fasta = ""
line_length = 70
genome_length_one_line = 0
CP_RLE = 0
CP_file_size = 0
genome_name = sys.argv[1][:-6]


genome_fasta_file = sys.argv[1]
min_k = int(sys.argv[2])
k = int(sys.argv[3])
partition = int(sys.argv[4])


n = 10
reads_length = 0
threshold = 1
genome_length = 0



def convert_genome_to_one_line(file, out_file):
    global genome_length, line_length, genome_length_one_line 

    f = open(out_file, 'w')
    fin = open(file , "r")

    i = -1

    for line in fin:
	i += 1
	if line[0] == ">":
		if i == 0:
			f.write(line)
			genome_length += len(line)
		else:
			f.write("\n" + line)
			genome_length += len(line)
	else:
		f.write(line.strip().upper())
		genome_length += len(line)
		genome_length_one_line += len(line.strip())

	if i == 1:
		line_length = len(line.strip())

    f.close() 

def run_length_encode(file, out_file):

	global genome_length, genome_length_one_line, RLE_numbers_dict, CP_RLE

	g = 0
	f = open(out_file, 'w')
	fin = open(file , "r")
	for line in fin:
		current = line[0]
        	counter = 1
		l = len(line)
		if line[0] == ">":
			f.write(line)
		else:
			for base in line[1:]+"$":
				g += 1
				if base == current:
                			counter += 1
        	    		else:
					if current == "A" or current == "C" or current == "G" or current == "T":		
						if counter > 10:
							f.write(current + str(counter))
							RLE_numbers_dict[str(counter)] += 1
						else:
							f.write(current*counter)				
					else:
						if counter > 2:
							f.write(current + str(counter))
							RLE_numbers_dict[str(counter)] += 1
						else:
							f.write(current*counter)


        	        		current = base
			                counter = 1
	
	f.close()

	# compute the genoe size after run-length encoding

	RLE_genome_length = 0
	fin = open(out_file, "r")
    	for line in fin:
        	for base in line:
			if base =="0" or base =="1":
				RLE_genome_length += 1
			else:
				RLE_genome_length += 8

	CP_RLE += RLE_genome_length/float(g *8)

			




def partition_genome_sequence(p, file):

	num_lines = sum(1 for line in open(file))
 
	size_part = num_lines/p
	if size_part < 2:
		print "\nError: Each partition must contain at least 2 lines"
		return 
	


	k = -1
	i = 1
	f = open(file[:-6]+"_"+str(i)+".fasta", 'w')
	fin = open(file , "r")
	for line in fin:
		k += 1
		f.write(line)
		if k%size_part == 0 and k != 0 :
			if i != p:
				f.close()
				i += 1
				f = open(file[:-6]+"_"+str(i)+".fasta", 'w')






def compute_freq_kmer(min_k, kmer_size, file):

    global kmer_dict, bases_freq_dict

    fin = open(file, 'r')
    for line in fin:
	l = len(line)
    	for i in range(l):
		for j in range(min_k, kmer_size + 1):
				wind = line[i: i + j]
				if i+j > l:
					if len(wind) == j:
						kmer_dict[wind] += 1
        			else:
					kmer_dict[wind] += 1

    #print "collecting kmers frequencies finished"



def compute_non_DNA_base_freq(file):

    fin = open(file, 'r')
    for line in fin:
	if line[0] != ">":
	    	for base in line.strip():
			if not base.isdigit():
				if base == 'A' or base == 'C' or base == 'T' or base == 'G':
					pass
				else:
					None_DNA_bases_freq_dict[base] += 1


    #print "collecting single bases frequencies finished"



# compress the genome using the codewords
def compress_reads(file1, out_file, codes_list, kmer_size):
    # decrease k by one becuase it will increased for the for loop

    global compressed_genome, genome_length, min_k, CP, line_length

    #print "start compressing"
    #classify thecode_words
    codeword_dict = defaultdict(str)

    for codeword in codes_list:
        codeword_dict[codeword[0]]=codeword[1]



    f_genome = open("temp.bits", 'w')


    fin = open(file1, 'r')
    for line in fin:
		if line[0] == ">":
			f_genome.write(line)

		else:			
			ll = len(line)
			i = 0
			while i < ll:
					comp_flag = False
					for j in range(kmer_size, 0, -1):
							if i + j <= ll:
								t = line[i:i+j]
								if t in codeword_dict:
									f_genome.write(codeword_dict[t])
									comp_flag = True
									i += j
									break


					if not comp_flag:
						tt = line[i]
						f_genome.write(tt)
						i += 1


    f_genome.close()

#convert the file of bits stings to ascii file and add the left over for each chromosome

    fin = open("temp.bits", 'r')
    f_genome = open(out_file, 'w')
    
    
    # add the codewords list as header
    # partition the code based on letters and integers to be easily parsed during uncompressing	

    RLE_codes = []
    huff_codes = []
    for i in range(len(codes_list)):
	if codes_list[i][0].isdigit():
		RLE_codes.append(codes_list[i])
	else:
		huff_codes.append(codes_list[i])
		
    header = ""
    for item in huff_codes:
	header += item[0]+item[1]

    # add ; to seperate the alphabitical and numerical codewords
    header += "+"

    
    if len(RLE_codes) > 0:
	for item in RLE_codes:
		header += item[0]+","+item[1] + ";"
        header = header[:-1]



    # add the priting line length
    header += ":" + str(line_length)
    f_genome.write("#UHT1_header:" + header)


    c = 0   # used to check of there is any header in the fasta file, if not create line that will store the leftover
    left_over = ""
    p = 0
    for line in fin:
	if line[0] == ">":
		left_over += "\n#UHT1_leftover:" + line.strip() + ":" + str(p)
		c += 1
	else:
		l = len(line)
		i = 0
		while i < l :
			if l - i < 8:
				left_over += ":" + line[i:].strip()
				break
			else:
				i += 8
				p += 8
    if c != 0:
	f_genome.write(left_over + "\n")

    if c == 0:
	f_genome.write("\n#UHT1_leftover:" + "None:-1:" + line[i:].strip() + "\n")

    fin = open("temp.bits", 'r')
    for line in fin:
	if line[0] != ">":
		l = len(line)
		i = 0
		while i < l :
			if l - i < 8:
				break
			else:
				f_genome.write(chr(int(line[i:i+8], 2)))
				i += 8

    f_genome.close()




# This function take as input dictionary where kmers are the IDs of the dictionary and the frequency associated with each kmer as the value of the dictionary
# and return anther dictionary where the frequency is the ID and the kmer/s that have the frequency as the value/s
def sort_kmers_dict_based_freq(dict):
	temp = defaultdict(list)

        for kmer in dict.keys():
            freq = dict[kmer]
            temp[freq].append(kmer)
        return temp


# This function take a dictionary where the frequency is the ID of the dictionary, and return the kmers with the top n frequencies
def get_n_max_freq_kmers(dict, n):
    a = []
    j = 0
    for i in sorted(dict.keys(), reverse=True):
		if j < n:
			a.append((i, dict[i][0]))
			#print "**", i, dict[i][0]
			j += 1
		else:
			pass
    return a


# This function for finding the kmers set that force the hufman tree to be unbalanced
def get_kmers_matching_hufman_design(dict):
    a = []

    keys = sorted(dict.keys())
    l = len(keys)
    if l <= 1:
	#print "no kmer freq collected"
	return a

    # compute the average
    sum = 0
    for e in keys:
        sum += int(e)

    temp_dict = dict


    keys = sorted(temp_dict.keys())


    new = sum
    r = 1

    l = len(keys)
    for i in range(l-1, 0, -1):
	t = temp_dict[keys[i]][0]

	if keys[i] < new/r :
		    f = True
		    for d in a:
			tt = len(d[1])
			if d[1] in t:
				f = False




		    if f:
            		a.append((keys[i], t))
	         	new = keys[i]
			r *= 2


    return a


# Huffman Coding #########################
# dont wory about these functions
def makeHuffTree(symbolTupleList):
   trees = list(symbolTupleList)

   heapq.heapify(trees)
   while len(trees) > 1:
      childR, childL = heapq.heappop(trees), heapq.heappop(trees)
      parent = (childL[0] + childR[0], childL, childR)
      heapq.heappush(trees, parent)
   return trees[0]

def create_Huff_code(huffTree, prefix = ''):
   global code_bits
   if len(huffTree) == 2:
      code_bits.append((huffTree[1], prefix))
   else:
      create_Huff_code(huffTree[1], prefix + '0')
      create_Huff_code(huffTree[2], prefix + '1')
   return code_bits




###### The main finction where we run the flow of the algorithm ########################

t1 = datetime.now()
num_kmer_to_encode = int(sys.argv[5])

partition_genome_sequence(partition, genome_fasta_file)

for i in range(1, partition+1):


	convert_genome_to_one_line(genome_name + '_'+str(i)+'.fasta', genome_name + '_one_line_'+str(i)+'.fasta')
	os.remove(genome_name + '_'+str(i)+'.fasta')

	run_length_encode(genome_name + '_one_line_'+str(i)+'.fasta', genome_name + '_one_line_RLE_'+str(i)+'.fasta')
	os.remove(genome_name + '_one_line_'+str(i)+'.fasta')

	compute_non_DNA_base_freq(genome_name + '_one_line_RLE_'+str(i)+'.fasta')
	Non_DNA_bases_freq = []
    	for key in sorted(None_DNA_bases_freq_dict.keys()):
		Non_DNA_bases_freq.append((None_DNA_bases_freq_dict[key], key))

	RLE_numbers_freq = []
    	for key in sorted(RLE_numbers_dict.keys()):
		RLE_numbers_freq.append((RLE_numbers_dict[key], key))

	compute_freq_kmer(min_k, k, genome_name + '_one_line_RLE_'+str(i)+'.fasta')
	kmer_sorted_freq_dict = sort_kmers_dict_based_freq(kmer_dict)

	# compression stage


	kmers_list = get_n_max_freq_kmers(kmer_sorted_freq_dict, num_kmer_to_encode) + RLE_numbers_freq + Non_DNA_bases_freq

	huffTree  = makeHuffTree(kmers_list)
	final_code_bits = create_Huff_code(huffTree)
	#print "--- Codewords used", final_code_bits

	compress_reads(genome_name + '_one_line_RLE_'+str(i)+'.fasta', genome_name + '_'+str(i)+'.uht', final_code_bits, k)
	CP_file_size += os.path.getsize(genome_name + '_'+str(i)+'.uht')

	code_bits = []
	Non_DNA_bases_freq = []

        RLE_numbers_dict.clear()
	kmer_dict.clear()
	kmer_sorted_freq_dict.clear()
	bases_freq_dict.clear()
	None_DNA_bases_freq_dict.clear()

	os.remove(genome_name + '_one_line_RLE_'+str(i)+'.fasta')


os.remove('temp.bits')



print "\nCompression ratio after run-length encoding:", CP_RLE/float(partition)

# compute file sizes before and after, Compression ratio

print "Final compression ratio based on file sizes, original and compressed:", CP_file_size/float(os.path.getsize(genome_fasta_file)) * 100

print "Compression time in seconds", (datetime.now() - t1).seconds, "\n"


