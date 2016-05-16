
from itertools import groupby
import string
from collections import defaultdict
from datetime import datetime
import sys
import os

compressed_file = sys.argv[1]
fil_name = sys.argv[1][:-4]


leftover_dict = defaultdict(list)
codes_list = []
line_length = 0


def convert_compressed_file_to_bits_file(file, out_file, marker):
    global codes_list, leftover_dict, line_length

    f = open(out_file, 'w')
    fin = open(file , "r")

    i = -1

    #collect leftover
    for line in fin:
	if "#UHT1_leftover" in line:
		tmp = line.strip().split(":")
		if len(tmp) == 4:
			leftover_dict[tmp[2]].append(tmp[1])
			leftover_dict[tmp[2]].append(tmp[3])
		if len(tmp) == 3:
			leftover_dict[tmp[2]].append(tmp[1])
			leftover_dict[tmp[2]].append("")



	
    fin.close() 
    
    fin = open(file , "r")
    pp = 0
    p = 0	
    last_key = ""
    for line in fin:
	i += 1
	if line[0:12] == "#UHT1_header":
		tmp = line.split(":")
		if marker == 1:
			line_length = int(tmp[2])

		# collect the code words
		a = [''.join(g) for _, g in groupby(tmp[1].split("+")[0], str.isalpha)]
		
		i = 0
		while i < len(a) - 1:
			codes_list.append((a[i], a[i+1]))
			i += 2	

		a = tmp[1].split("+")[1].split(";")	
		if len(a) > 0:
			for item in a:
				t = item.split(",")
				if len(t) > 1:
					codes_list.append((t[0], t[1]))
		
	elif line[0:14] != "#UHT1_leftover":
		for c in range(len(line)):			 
			if str(p) in leftover_dict:
				if pp == 0 and "-1" in leftover_dict:
					f.write(leftover_dict['-1'][1] + "\n")
					f.write(leftover_dict[str(p)][0] + "\n")
 
					b = bin(ord(line[c]))[2:].zfill(8) 
					f.write(b)
					last_key = str(p) 

				elif pp == 0:
					f.write(leftover_dict[str(p)][0] + "\n")

					b = bin(ord(line[c]))[2:].zfill(8) 
					f.write(b)
					last_key = str(p) 
				else:
					f.write(leftover_dict[last_key][1] + "\n")
					f.write(leftover_dict[str(p)][0] + "\n")

					b = bin(ord(line[c]))[2:].zfill(8) 
					f.write(b)
					last_key = str(p) 
				pp += 1

			else:
				b = bin(ord(line[c]))[2:].zfill(8) 
				f.write(b)
			p += 8

			
    if len(leftover_dict.keys()) == 1 and "-1" in leftover_dict:
    	f.write(leftover_dict["-1"][1])
    elif len(leftover_dict.keys()) == 1 and last_key == "" and "-1" not in leftover_dict:
	for key in leftover_dict.keys():
		f.write(leftover_dict[key][0] + "\n")
		f.write(leftover_dict[key][1] )
    else:
	f.write(leftover_dict[last_key][1])
    
	

    f.close() 



# uncompress the compressed genome using the codes list
def decompress_genome_bits(file, out_file, codes_list, part_num):
    global line_length

    #print "\nstart decompressing"

    bitword_dict = defaultdict(str)
    for codeword in codes_list:
        bitword_dict[codeword[1]]=codeword[0]

    #print codes_list

    # find the maximum length of bitword
    max_bits = 0
    min_bits = 1000
    for codeword in codes_list:
        length = len(codeword[1])
        if length > max_bits:
            max_bits = length
        if length < min_bits:
            min_bits = length

    last_base = ""
    temp_print_line = ""
    f = open(out_file, 'a')
    fin = open(file , "r")
    pos = 0
    for line in fin:
        i = 0
	l = len(line)
	if line[0] == ">":
		if part_num == 1:
			if pos == 0:
				f.write(line)		
				pos += 1
			else:
				f.write("\n" + line)		
		else:
			f.write("\n" + line)
			pos += 1		
	
	else:
	    if pos == 0 and part_num > 1:
		f.write("\n")

	    temp_print_line = ""	
	    while i < l:	
	            comp_flag = False
        	    for j in range(max_bits, min_bits - 1 , -1):
                	if i + j <= l:
	                    t = line[i : i + j]
	                    if t in bitword_dict:
				if bitword_dict[t].isdigit():
				    temp_print_line += (int(bitword_dict[t]) - 1) * last_base
                	            comp_flag = True
				    i += j
	                            break
				else:
				    temp_print_line += str(bitword_dict[t])
                	            comp_flag = True
				    i += j
	                            break
         			
            	    if not comp_flag:
			temp_print_line += line[i]
			i += 1

 		    last_base = temp_print_line[-1]

		    while len(temp_print_line) > line_length:
			f.write(temp_print_line[:line_length] + "\n")
			temp_print_line = temp_print_line[line_length:]	

	    f.write(temp_print_line.strip())
    


    f.close()



#def decompress_genome_RLE(comp_genome, codes_list):

t1 = datetime.now()

i = 1
flag = True

# remove old file so no appending occure on the old files
if os.path.isfile(fil_name + ".fasta"):
	os.remove(fil_name + ".fasta")

while flag:

	if str(i) == "1":
		convert_compressed_file_to_bits_file(compressed_file, fil_name + ".bits", 1)
		os.remove(compressed_file)

		decompress_genome_bits(fil_name + ".bits", fil_name + ".fasta", codes_list, 1)
		os.remove(fil_name + ".bits")
		
	elif os.path.isfile(fil_name[:-1] + str(i) + ".uht"):
		compressed_file = fil_name[:-1] + str(i) + ".uht"

		convert_compressed_file_to_bits_file(compressed_file, compressed_file[:-5] + str(i) + ".bits", i)
		os.remove(compressed_file)

		decompress_genome_bits(compressed_file[:-5] + str(i) + ".bits", fil_name + ".fasta", codes_list, i)
		os.remove(compressed_file[:-5] + str(i) + ".bits")

	else:
		flag = False
	
	
	i += 1
	leftover_dict.clear()
	codes_list = []
			


print "Decompression has finished in ", (datetime.now() - t1).seconds, "seconds\n"














