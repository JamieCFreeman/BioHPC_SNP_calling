
#!/bin/bash

# JCF May2019


# Must be run in parent directory that has subfolders: 
# 		dedup_bam, bam & bai files
#		gvcf, where g.vcf files .tbi index & .log file will output	

# Requires variables:
# REF="/workdir/jcf236/bwa_dir/mds_ref.fa.gz", reference fasta as .fa !!!! GATK won't take .gz
# CORES to set # threads, bwa uses as many as you give it
# ID is sample name prefix for file
# INT is interval "chr1" seems to work
# List of scaffolds obtained from the fai file
	# awk '{print $1}' mds_ref.fa.fai  > /workdir/jcf236/Workflow_biohpc/IL_scaff_list.txt

# On BioHPC: Version:4.0.1.1
# On my desktop:
# export PATH=$PATH:/home/Jamie/data/Genomics_programs/gatk-4.1.1.0


# There are pipes in the Illunima assembly scaffold names- can't go into file names
INT_short=$( echo "$INT" | sed -n -e 's/.*ref|//p' | sed -n 's/|$//p' )



/programs/gatk4/gatk --java-options "-Xmx4g" HaplotypeCaller  \
-R ${REF} \
-I ./dedup_bam/${ID}_dedup.bam \
-O ./gvcf/${ID}_${INT_short}.g.vcf.gz \
-L ${INT} \
-ERC GVCF 2> ./gvcf/${ID}_${INT_short}.log
   




