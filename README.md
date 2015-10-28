
HGA tool
This tool helps to apply the Hierarchical Genome Assembly (HGA) method.

--Software version--
This is version 1.0.0, supporting complete HGA assembly using SPAdes and velvet to assemble the partitions, for now.

--Prerequisite--
Velvet should be installed with option [make 'MAXKMERLENGTH=111' 'LONGSEQUENCES=1']. 'LONGSEQUENCES=1' to allow Velvet to accept contigs (long sequences) as input; while 
'MAXKMERLENGTH=111' to setup assembly using kmer of size up tp 111 as the default max kmer length is 31.

--Installation--
The tool is built using python; there is no need for installation.

--Running HGA Tool--
Run the command 
python HGA.py [options]

where the options are:
--Parameters--:
 -velvet  path      Path to the directory where velvet bineries (velveth, velvetg) are located
 -spades  path      Path to the directory where spades.py file is located
 -PA      string    select one of {SPAdes, velvet} that will assemble the partitions
 -P12     file      interlaced fastq file of paired reads that will be used in the partiton step, usually raw or clean (corrected) reads
 -R12     file      interlaced fastq file of paired reads that will be used in the re-assembly step, usually raw or clean (corrected) reads
 -ins     int       insert size of the fragments
 -std     int       standard deviation of the inset size
 -P       int       number of partitions
 -Pkmer   int(odd)  kmer value for assembling the partitions
 -Rkmer   int(odd)  kmer value for re-assembly step
 -t       int       number of threads to be used for re-assembly step using SPAdes assembler default 1
 -out     Path      output path
 -h                 print option menue

-Output formats
After HGA finish running, the following will be at the given output directory:
	- All, 1-p, partitioned fastq files; and their assemblies each in folder under the name of the partition.
	- ./combined_contigs/combined_contigs.fasta; fasta file of merging all contigs of the partitions assemblies.
	- ./merged_contigs/merged_contigs.fasta, fasta file of combining all contigs of the partitions assemblies.
	- ./HGA_merged/contigs.fasta (as well scaffolds.fasta); contigs of HGA(merged contigs) flow.
	- ./HGA_combined/contigs.fasta (as well scaffolds.fasta); contigs of HGA(combined contigs) flow.


--Contacts--
For more quiries or questions, please contact
Anas Al-okaily, aaa10013@engr.uconn.edu


Last update: 10/2015
