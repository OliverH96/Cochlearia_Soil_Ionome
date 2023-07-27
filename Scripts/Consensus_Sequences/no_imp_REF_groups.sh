#!/bin/bash
#SBATCH --job-name=topSNP_REF
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24g
#SBATCH --time=024:00:00
#SBATCH --output=/gpfs01/home/mbxoh2/OandE/%x.out
#SBATCH --error=/gpfs01/home/mbxoh2/OandE/%x.err

source $HOME/.bash_profile
cd ~/genome_data/cochlearia/HKT1_consensus
conda activate gatk

ref=~/genome_data/reference/C_excelsa_V5.fasta
refvcf=no_imp_hkt1_REF.vcf.gz

##Top rated SNP REF
#Group1 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn BNK_21 -sn CHT_2 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group1_REF.fa
rm top-rated-overall_REF.vcf*

#Group2 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn BRE_1 -sn BRE_2 -sn BRE_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group2_REF.fa
rm top-rated-overall_REF.vcf*

#Group3 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn CHA_1 -sn CHA_2 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group3_REF.fa
rm top-rated-overall_REF.vcf*

#Group4 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group4_REF.fa
rm top-rated-overall_REF.vcf*

#Group5 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn KVA_10 -sn KVA_2 -sn KVA_3 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group5_REF.fa
rm top-rated-overall_REF.vcf*

#Group6 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn LAB_1 -sn LAB_2 -sn LAB_3 -sn LAB_4 -sn LAB_5 -sn LAB_6 -sn LAB_7 -sn LAB_8 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group6_REF.fa
rm top-rated-overall_REF.vcf*

#Group7 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn LNL_1 -sn LNL_8 -sn OSW_1 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group7_REF.fa
rm top-rated-overall_REF.vcf*

#Group8 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn NEN_1 -sn NEN_200 -sn NEN_3 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group8_REF.fa
rm top-rated-overall_REF.vcf*

#Group9 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group9_REF.fa
rm top-rated-overall_REF.vcf*

#Group10 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn SAL_3 -sn SAL_4 -sn SAL_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group10_REF.fa
rm top-rated-overall_REF.vcf*

#Group11 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $refvcf --select "AF > 0.5" -sn TRO_1 -sn TRO_5 -sn TRO_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_REF.vcf.gz -o  top_rated_SNP/top-rated-overall-group11_REF.fa
rm top-rated-overall_REF.vcf*
