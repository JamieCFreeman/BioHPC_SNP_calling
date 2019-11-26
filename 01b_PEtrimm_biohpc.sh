
#!/bin/bash

# 2019-11-26 JCF

# General command for using Trimmomatic to trim single end reads for quality & adaptor sequences
# Must be run in parent directory that has subfolders: 
#		fastq, for file to trim
# 		trim_fastq, for output
	

# Requires variables:
# CORES to set # threads, trimmomatic seems to use ~1.5 no matter assignment
# FILE is input file name (output from sequencing center has long string, usually want to rename to something shorter)
# ID is sample name prefix for output file

# Trimmomatic version trimmomatic-0.39 (checked 2019-5-22) 
# Parameters from http://www.usadellab.org/cms/?page=trimmomatic PE (2019-11-26)

# Trimmomatic doesn't ever seem to use >2 cores

java -jar /programs/trimmomatic/trimmomatic-0.39.jar PE -threads ${CORES} -phred33 \
	./fastq/${ID}_1.fastq ./fastq/${ID}_2.fastq \
	./trim_fastq/${ID}_trim_P_1.fastq.gz ./trim_fastq/${ID}_trim_UN_1.fastq.gz \
	./trim_fastq/${ID}_trim_P_2.fastq.gz ./trim_fastq/${ID}_trim_UN_2.fastq.gz \
	ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:36 2> ./trim_fastq/${ID}_trim.log

# Set exit status as status of last command
exit $?
