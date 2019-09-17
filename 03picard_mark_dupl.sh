#!/bin/bash

# JCF May2019

# General command for using Picard to mark duplicate reads in a bam file.
# Syntax appropriate for Picard 2.19.2, they will change syntax some time next year 

# Must be run in parent directory that has subfolders: 
#		sorted_bam, the source of the files to be marked
#		dedup_bam, the directory where the output bam & index will be writen

# Requires variables:
# ID is sample name prefix for file

# BioHPC Picard jar path: /programs/picard-tools-2.19.2/picard.jar
# Desktop Picard jar path: /home/jamie/Genomics_programs/picard/picard-2.19.2.jar 
# Also have: /home/jamie/Genomics_programs/picard/picard-2.20.1.jar (Latest version as of 2019-5-22)


java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates \
	INPUT="./sort_bam/${ID}_sort.bam" \
	OUTPUT="./dedup_bam/${ID}_dedup.bam" \
	METRICS_FILE="./dedup_bam/${ID}_sorted_dedup_metrics.txt"
	samtools index "./dedup_bam/${ID}_dedup.bam" 2> "./dedup_bam/${ID}_markdups.log"



# 2> doesn't deposit anything in log, unless command fails or warnings exist.
# Test whether log file is empty & delete if it is. 
# The expression -s is TRUE if the file exists and has a size > 0.
if [ -s "./dedup_bam/${ID}_markdups.log" ] 
	then echo ""
	else rm "./dedup_bam/${ID}_markdups.log"
fi






