#!/bin/bash

# JCF May2019


# Must be run in parent directory that has subfolders: 
# 		dedup_bam, bam & bai files
#		gvcf, where g.vcf files .tbi index & .log file will output

# Requires variables:
# REF="/workdir/jcf236/bwa_dir/mds_ref.fa.gz", reference fasta as .fa !!!! GATK won't take .gz
# CORES to set # threads, bwa uses as many as you give it
# ID is sample name prefix for file
# MAP_PATH is path to txt file sample map of vcfs to be imported

# On BioHPC: Version:4.0.1.1
# On my desktop:
# export PATH=$PATH:/home/Jamie/Genomics_programs/gatk-4.1.1.0
#  awk '{print $1}' /workdir/jcf236/bwa_dir/mds_ref.fa.fai > /workdir/jcf236/bwa_dir/mds_scaff_list 

# In a loop, call variants for INT between num_start & num_end

# Make a temporary directory to store individual vcf files
mkdir "./vcf/temp${num_start}-${num_end}"

while read -r INT; do

    echo "$INT"
    # There are pipes in the Illunima assembly scaffold names- can't go into file names
    INT_short=$( echo "$INT" | sed -n -e 's/.*ref|//p' | sed -n 's/|$//p' )


    /programs/gatk4/gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      --sample-name-map ${MAP_PATH} \
      --genomicsdb-workspace-path ${GDB_PATH}/GDB_${INT_short} \
      --intervals ${INT} 2> "${GDB_PATH}/${INT_short}.log"

    /programs/gatk4/gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R ${REF} \
      -V gendb:///${GDB_PATH}/GDB_${INT_short} \
      -stand-call-conf 30 \
      -O "./vcf/temp${num_start}-${num_end}/jointcalls_${INT_short}.vcf" 2> "./vcf/temp${num_start}-${num_end}/jointcalls_${INT_short}.log"

    # Done with 
    rm -r ${GDB_PATH}/GDB_${INT_short}
    rm "${GDB_PATH}/${INT_short}.log"
    
    # Remove GDB_log files if the vcf exists
    #if [ -e "jointcalls_${num_start}-${num_end}.vcf"]
    #then
    #rm  "${GDB_PATH}/${INT_short}.log"
    #fi


done < <(printf '%s\n' "$INT_list")

cd ./vcf/temp${num_start}-${num_end}
# Get a list of the vcf files to be combined
find ~+ -maxdepth 1 -iname "jointcalls_*" -iname "*.vcf" | sort -V > "vcfs_to_combine.txt"

cd ..
# Use GatherVcfs to combine them (must be properly sorted!)
java -jar /programs/picard-tools-2.19.2/picard.jar GatherVcfs \
    I="./temp${num_start}-${num_end}/vcfs_to_combine.txt" \
    O="jointcalls_${num_start}-${num_end}.vcf"

    
    
