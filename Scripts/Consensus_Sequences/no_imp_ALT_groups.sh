#!/bin/bash
#SBATCH --job-name=topSNP_ALT
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
altvcf=no_imp_hkt1_ALT.vcf.gz

##Top rated SNP ALT
#Group1 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group1_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group2 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn ALO_17 -sn ALO_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group2_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group3 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn BEA_2 -sn BEA_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group3_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group4 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn BRI_2 -sn BRI_5 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group4_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group5 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn ELI_1 -sn ELI_2 -sn ELI_3 -sn ELI_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group5_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group6 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn FDE_2 -sn FDE_3 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group6_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group7 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn FLA_2 -sn FLA_3 -sn FLA_4 -sn FLA_5 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group7_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group8 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn FRE_13 -sn RYE_1 -sn SCO_1 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group8_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group9 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn FTW_1 -sn FTW_3 -sn FTW_5 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group9_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group10 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn GEO_2 -sn GEO_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group10_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group11 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn GRE_1 -sn GRE_2 -sn GRE_3 -sn GRE_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group11_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group12 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group12_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group13 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn LOS_1 -sn LOS_2 -sn LOS_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group13_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group14 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group14_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group15 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn ODN_1 -sn BRE_3 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group15_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group16 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn OSW_2 -sn OSW_3B -sn OSW_6 -sn OSW_7 -sn OSW_8 -sn OSW_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group16_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group17 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn ROT_13 -sn ROT_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group17_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group18 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group18_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group19 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn SKF_3 -sn SKF_5 -sn SKF_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group19_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group20 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn SKI_5 -sn VAG_7 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group20_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group21 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn SKN_1 -sn SKN_2 -sn SKN_5 -sn SKN_8 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group21_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group22 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn TET_2 -sn TET_4 -sn TET_6 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group22_ALT.fa
rm top-rated-overall_ALT.vcf*

#Group23 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $altvcf --select "AF > 0.5" -sn UIG_1 -sn UIG_10 -sn UIG_2 -sn UIG_3 -sn UIG_4 -sn UIG_5 -sn UIG_6 -sn UIG_7 -sn UIG_8 -sn UIG_9 -L Cexcelsa_scaf_2:4113454-4117809 -O top-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus top-rated-overall_ALT.vcf.gz -o  top_rated_SNP/top-rated-overall-group23_ALT.fa
rm top-rated-overall_ALT.vcf*
