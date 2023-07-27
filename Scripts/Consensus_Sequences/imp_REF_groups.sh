#!/bin/bash
#SBATCH --job-name=imptopSNP_REF
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
vcf=imp_hkt1.vcf.gz

##Top rated SNP REF
#Group1 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn AUV_1 -sn AUV_2 -sn AUV_3 -sn AUV_4 -sn AUV_5 -sn AUV_6 -sn AUV_7 -sn AUV_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group1_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group2 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn BEI_1 -sn BEI_10 -sn BEI_2 -sn BEI_3 -sn BEI_4 -sn BEI_5 -sn BEI_6 -sn BEI_7 -sn BEI_8 -sn BEI_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group2_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group3 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn BNK_21 -sn BRE_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group3_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group4 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn CHA_1 -sn CHA_2 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group4_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group5 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn DIE_1 -sn DIE_2 -sn DIE_3 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group5_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group6 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn DOT_1 -sn DOT_2 -sn DOT_4 -sn DOT_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o imp_top_rated_SNP/top-rated-group6_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group7 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GAL_1 -sn GAL_2 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group7_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group8 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GAN_1 -sn GAN_3 -sn GAN_4 -sn GAN_5 -sn GAN_6 -sn GAN_7 -sn GAN_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group8_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group9 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GAR_3 -sn GAR_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group9_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group10 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GOR_4 -sn GOR_5 -sn GOR_6 -sn GOR_7 -sn GOR_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group10_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group11 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn JAZ_1 -sn JAZ_2 -sn JAZ_3 -sn JAZ_4 -sn JAZ_5 -sn JAZ_6 -sn JAZ_7 -sn JAZ_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group11_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group12 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn JOR_12 -sn JOR_13 -sn JOR_3 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group12_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group13 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LAB_1 -sn LAB_2 -sn LAB_3 -sn LAB_4 -sn LAB_5 -sn LAB_6 -sn LAB_7 -sn LAB_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group13_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group14 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LAE_2 -sn LAE_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group14_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group15 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LAW_1 -sn LAW_3 -sn LAW_4 -sn LAW_5 -sn LAW_7 -sn LAW_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group15_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group16 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LNL_8 -sn LOC_1 -sn KAU_1 -sn DYR_1 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group16_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group17 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LUT_1 -sn LUT_2  -sn LUT_3 -sn LUT_4 -sn LUT_5 -sn LUT_6a -sn LUT_6b -sn LUT_7 -sn LUT_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group17_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group18 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn NEN_1 -sn NEN_200 -sn NEN_3 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group18_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group19 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group19_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group20 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn OSW_1 -sn RUN_7 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group20_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group21 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn PUY_1 -sn PUY_2 -sn PUY_3 -sn PUY_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group21_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group22 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn ROS_1 -sn ROS_2 -sn ROS_3 -sn ROS_4 -sn ROS_5 -sn ROS_6 -sn ROS_7 -sn ROS_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group22_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group23 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn SAL_3 -sn SAL_6 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group23_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group24 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn STE_1 -sn STE_2 -sn STE_3 -sn STE_5 -sn STE_6 -sn STE_7 -sn STE_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group24_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group25 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn VIL_7 -sn RUZ_1 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group25_REF.fa
rm imptop-rated-overall_REF.vcf*

#Group26 - filter vcf and generate REF sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn WOL_10 -sn WOL_2 -sn WOL_6 -sn WOL_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_REF.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_REF.vcf.gz -o  imp_top_rated_SNP/top-rated-group26_REF.fa
rm imptop-rated-overall_REF.vcf*
