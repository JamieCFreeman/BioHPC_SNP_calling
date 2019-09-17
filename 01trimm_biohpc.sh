#!/bin/bash

# 2019-5-22 JCF

# General command for using Trimmomatic to trim single end reads for quality & adaptor sequences
# Must be run in parent directory that has subfolders: 
#		fastq, for file to trim
# 		trim_fastq, for output
	

# Requires variables:
# CORES to set # threads, trimmomatic seems to use ~1.5 no matter assignment
# FILE is input file name (output from sequencing center has long string, usually want to rename to something shorter)
# ID is sample name prefix for output file

# Trimmomatic version trimmomatic-0.39 (2019-5-22) 

java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -threads ${CORES} -phred33 ./fastq/${FILE} \
	./trim_fastq/${ID}_trim.fastq.gz ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
	SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 2> ./trim_fastq/${ID}_trim.log

# Set exit status as status of last command
exit $?
