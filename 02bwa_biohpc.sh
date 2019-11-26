#!/bin/bash

# edit 2019-11-26

# General command for using bwa mem to map short reads to a reference genome
# Must be run in parent directory that has subfolders: 
# 		trim_fastq, including fastq.gz files for mapping
# 		bam, the destination folder for bam files, their index files, & their flagstat output
#		sorted_bam, the destination for the mapped & sorted bam files, their index files & their flagstat output	


# bwa version 0.7.13 as installed on BioHPC

# Requires variables:
# RG_string is read group information, ex "@RG\tID:BSAF7_1_RES_b\tSM=:SAF7_1_RES\tPU:none\tLB:D01\tPL:ILLUMINA"
# REF="/workdir/jcf236/bwa_dir/mds_ref.fa.gz", reference fasta as .fa or .fa.gz
# CORES to set # threads, bwa uses as many as you give it
# ID is sample name prefix for file
# PE is yes/no: yes when data is pair end, no when single 

# Check that variables are set & exit otherwise
if [ -z "$RG_string" ] || [ -z "$REF" ] || [ -z "$CORES" ] || [ -z "$ID" ] || [ -z "$PE" ];
then
 echo "One or more variables are undefined."
 exit 1
fi

# Single-end command:
if [ "$PE" = "no" ]; 
	then
	# Pipe bwa mem output directly to samtools to make bam (sam bigger)
	bwa mem -M -R ${RG_string} -t ${CORES} ${REF} ./trim_fastq/${ID}_trim.fastq.gz 2> ./bam/${ID}_bwa.log | \
		samtools view -bS - > ./bam/${ID}.bam
		samtools flagstat ./bam/${ID}.bam > ./bam/${ID}.stats
		samtools sort -@ ${CORES} ./bam/${ID}.bam > ./sort_bam/${ID}_sort.bam
		samtools index ./sort_bam/${ID}_sort.bam
		samtools flagstat ./sort_bam/${ID}_sort.bam > ./sort_bam/${ID}_sort.stats
	fi

# Paired-end command
if [ "$PE" = "yes" ]; 
	then
	# Pipe bwa mem output directly to samtools to make bam (sam bigger)
	bwa mem -M -R ${RG_string} -t ${CORES} ${REF} ./trim_fastq/${ID}_trim_P_1.fastq.gz \
		./trim_fastq/${ID}_trim_P_2.fastq.gz  2> ./bam/${ID}_bwa.log | \
		samtools view -bS - > ./bam/${ID}.bam
		samtools flagstat ./bam/${ID}.bam > ./bam/${ID}.stats
		samtools sort -@ ${CORES} ./bam/${ID}.bam > ./sort_bam/${ID}_sort.bam
		samtools index ./sort_bam/${ID}_sort.bam
		samtools flagstat ./sort_bam/${ID}_sort.bam > ./sort_bam/${ID}_sort.stats
	fi

# Compare sorted stats file to unsorted stats file. Should be the same- # of reads will only change if one is truncated.
# If they are the same, we want to delete the unsorted.
# If they aren't the same, keep both files
cmp -s ./bam/${ID}.stats ./sort_bam/${ID}_sort.stats
if cmp -s ./bam/${ID}.stats ./sort_bam/${ID}_sort.stats; then
	rm ./bam/${ID}.bam;
fi
