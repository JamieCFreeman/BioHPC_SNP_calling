#!/bin/bash

# 2019-5-22 JCF

# Must exist in analysis directory:
#	./fastq folder with .fq files to be analyzed
#	Folder with reference genome & associated index files for bwa mem & GATK. GATK needs unzipped .fa.
#			Example commands for creating index files below.
# 	The reference genome folder should also contain a file that contains the list of scaffold names in
#			the reference genome. This is used to scatter jobs by scaffold at the gvcf & vcf steps.
#			Example commands to create file below.  Extra caution is required when using NCBI ref genomes,
# 			as the scaffold names can have illegal characters for file names.

# Should be folder containing ./fastq & also the reference assembly directory
# Remember that GATK needs uncompressed reference file 


# 3 levels of setup
# 3 (whole script): creates all index files, runs fastQC, creates directory structure for analysis, logs program version info
# 2 creates directory structure for analysis, logs program version info
# 1 Just logs program version info


# Log program version info

mkdir -p ./metadata

OUTPUT_FILE=$( echo  $( date +%Y-%m-%d )_$( echo $PWD | sed 's /.*/  g' )_version.log ) 

echo $( date ) > ./metadata/$OUTPUT_FILE
echo $PWD >> ./metadata/$OUTPUT_FILE

/programs/parallel/bin/parallel --version >> ./metadata/$OUTPUT_FILE
/programs/FastQC-0.11.8/fastqc --version >> ./metadata/$OUTPUT_FILE
samtools --version >> ./metadata/$OUTPUT_FILE
java -jar /programs/picard-tools-2.19.2/picard.jar CreateSequenceDictionary --version >> ./metadata/$OUTPUT_FILE
java -jar /programs/trimmomatic/trimmomatic-0.39.jar -version >> ./metadata/$OUTPUT_FILE
cmp --version >> ./metadata/$OUTPUT_FILE
/programs/gatk4/gatk HaplotypeCaller --version >> ./metadata/$OUTPUT_FILE

# If level = 1, then exit here.
if [ "$1" -eq 1 ]; then
 exit 0;
fi


# Set up environment with child directories for results. 
mkdir ./trim_fastq
mkdir ./bam
mkdir ./sort_bam
mkdir ./dedup_bam
mkdir ./gvcf_int
mkdir ./gvcf
mkdir ./vcf
mkdir ./GDB_import

# If level = 2, then exit here.
if [ "$1" -eq 2 ]; then
 exit 0;
fi


# Do you want to run fastqc? Will check quality metrics for fq data.

cd fastq
#find -iname "*.fastq.gz"

# Use GNU parallel to run all fastqc jobs
export PATH=/programs/parallel/bin:$PATH
find . -maxdepth 1 -iname "*.fastq.gz" | parallel -n 1 \
	/programs/FastQC-0.11.8/fastqc {}


mkdir ./fastqc
mv *.zip ./fastqc
mv *.html ./fastqc




# Do you need index files?
# Set variable for reference genome name/path
REF=round5.curated.fasta

# summary files needed for GATK, only need to be done once for reference sequence
samtools faidx ${REF}


bwa index -a bwtsw ${REF}
java -jar /home/jamie/Genomics_programs/picard.jar CreateSequenceDictionary \
R=/home/jamie/data/"bwa dir"/${REF} \
O=/home/jamie/data/"bwa dir"/mds_ref.dict


java -jar /programs/picard-tools-2.19.2/picard.jar CreateSequenceDictionary \
R=${REF} \
O=curated.dict



# When jobs are scattered (gvcf & vcf steps), a list of scaffold names is required
awk '{print $1}' ${REF}.fai > ./round5_ref_dir/round5.curated.scaff.list.txt

# Ref genome from NCBI has pipes in the scaffold names, need to pull them out.
awk '{print $1}' ${REF}.fai  | sed -n -e 's/.*ref|//p' | sed -n 's/|$//p' > /workdir/jcf236/mapped2IL/bwa_dir/mds_ref.scaff.list.txt








