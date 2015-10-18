
HGA tool
This tool helps to apply the Hierarchical Genome Assembly (HGA) method.

--Software version--
This is version 1.0.0, supporting complete HGA assembly using SPAdes and velvet to assemble the partitions, for now.

--Prerequisite--
Velvet should be installed with option 'LONGSEQUENCES=1' in the make command, to allow Velvet to accept contigs (long sequences) as input; as well and 
'MAXKMERLENGTH=111' as the default max kmer length is 31; and 'OPENMP=1' to allow the use of multiple CPU cores on the same machine.

--Installation--
The tool is built using python; there is no need for installation.

--Running HGA Tool--
Run the command 
python HGA.py [options]

where the options are:
--Parameters--:
 -velvet  path      Path to the velvet bineries
 -spades  path      Path to the directory where spades.py file
 -PA      string    select one of {SPAdes, velvet} that will assemble the partitions
 -12     file      fastq file of the interlaced paired reads
 -ins     int       insert size of the fragments
 -std     int       standard deviation of the inset size
 -P       int       number of partitions
 -Pkmer   int(odd)  kmer for assembling the partions
 -Rkmer   int(odd)  kmer value for re-assembly step
 -t       int       number of threads to be used for re-assembly step using SPAdes assembler
 -out     Path      output path
 -h                 print option menue

-Output formats
After HGA finish running, the following will be at the given output directory:
	- All, 1-p, partitioned fastq files.
	- ./combined_contigs/combined_contigs.fasta; fasta file of merging all contigs of the partitions assemblies.
	- ./merged_contigs/merged_contigs.fasta, fasta file of combining all contigs of the partitions assemblies.
	- ./HGA_merged/contigs.fasta (as well scaffolds.fasta); contigs of HGA(merged contigs) flow.
	- ./HGA_combined/contigs.fasta (as well scaffolds.fasta); contigs of HGA(combined contigs) flow.


--Contacts--
For more quiries or questions, please contact
Anas Al-okaily, aaa10013@engr.uconn.edu


Last update: 10/2015
