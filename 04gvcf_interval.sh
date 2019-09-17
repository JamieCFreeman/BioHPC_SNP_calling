#!/bin/bash

# JCF May2019

# Must be run in parent directory that has subfolders: 
# 		dedup_bam, bam & bai files
#		gvcf_int, where g.vcf files .tbi index & .log file will output	

# Requires variables:
# REF="/workdir/jcf236/bwa_dir/mds_ref.fa.gz", reference fasta as .fa !!!! GATK won't take .gz
# CORES to set # threads, bwa uses as many as you give it
# ID is sample name prefix for file
# INT is interval

# On BioHPC: Version: 4.0.1.1
# On my desktop:
# export PATH=$PATH:/home/Jamie/data/Genomics_programs/gatk-4.1.1.0


/programs/gatk4/gatk --java-options "-Xmx4g" HaplotypeCaller  \
-R ${REF} \
-I ./dedup_bam/${ID}_dedup.bam \
-O ./gvcf_int/${ID}_${INT}.g.vcf.gz \
-L ${INT} \
-ERC GVCF 2> ./gvcf_int/${ID}_${INT}.log
   
# Set exit status as status of last command
exit $?
