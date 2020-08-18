#!/bin/bash

# 2020-3-16 JCF

# Separate SNPs and indels from raw vcf file, filter them separately, and output vcf of PASS.

# Input:
#	./vcf/{VCF_pref}.vcf
#	${REF}

# Output:
#       ./vcf/{VCF_pref}_SNPs.vcf
#       ./vcf/{VCF_pref}_SNPs_filtered.vcf	
#       ./vcf/{VCF_pref}_SNPs_PASS.vcf
#       ./vcf/{VCF_pref}_INDELs.vcf
#       ./vcf/{VCF_pref}_INDELs_filered.vcf
#       ./vcf/{VCF_pref}_INDELs_PASS.vcf

# Filter SNPs
/programs/gatk4/gatk --java-options "-Xmx4g" SelectVariants \
	-R ${REF} \
	-V ./vcf/${VCF_pref}.vcf \
	--select-type-to-include SNP \
	-O ./vcf/${VCF_pref}_SNPs.vcf

/programs/gatk4/gatk --java-options "-Xmx4g"  VariantFiltration \
	-R ${REF} \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	--filter-name "Low_conf" -cluster 3 -window 10 \
	-V "./vcf/${VCF_pref}_SNPs.vcf" \
	-O "./vcf/${VCF_pref}_SNPs_filtered.vcf"

/programs/gatk4/gatk --java-options "-Xmx4g" SelectVariants \
	-R ${REF} \
	-V "./vcf/${VCF_pref}_SNPs_filtered.vcf" \
	--exclude-filtered \
	-O "./vcf/${VCF_pref}_SNPs_PASS.vcf"

# Filter indels


/programs/gatk4/gatk --java-options "-Xmx4g" SelectVariants \
	-R ${REF} \
	-V "./vcf/${VCF_pref}.vcf" \
	--select-type-to-include INDEL \
	-O "./vcf/${VCF_pref}_INDELs.vcf" 

	# Decimal points are necesaary to define as float instead of integer!
	#

/programs/gatk4/gatk --java-options "-Xmx4g"  VariantFiltration \
	-R ${REF} \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-filter "AC < 20.0" --filter-name "AC20" \
	-V "./vcf/${VCF_pref}_INDELs.vcf" \
	-O "./vcf/${VCF_pref}_INDELs_filtered.vcf" 

/programs/gatk4/gatk --java-options "-Xmx4g" SelectVariants \
	-R ${REF} \
	-V "./vcf/${VCF_pref}_INDELs_filtered.vcf" \
	--exclude-filtered \
	-O "./vcf/${VCF_pref}_INDELs_PASS.vcf"
	 
