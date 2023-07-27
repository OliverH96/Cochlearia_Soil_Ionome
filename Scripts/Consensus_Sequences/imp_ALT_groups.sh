#!/bin/bash
#SBATCH --job-name=imptopSNP_ALT
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

##Top rated SNP ALT
#Group1 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group1_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group2 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn ALO_17 -sn ALO_6 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group2_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group3 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn BRE_3 -sn ODN_1 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group3_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group4 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn BRI_2 -sn BRI_5 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group4_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group5 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn CHT_3 -sn CHT_5 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group5_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group6 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn CUM_1 -sn FOR_1 -sn FRE_13 -sn JON_1 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group6_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group7 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn DAR_1 -sn DAR_3 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group7_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group8 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn ELI_1 -sn ELI_2 -sn ELI_3 -sn ELI_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group8_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group9 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn FDE_2 -sn FDE_3 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group9_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group10 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn FLT_1 -sn FLT_3 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group10_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group11 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn FTW_1 -sn FTW_3 -sn FTW_5 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group11_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group12 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GAL_4 -sn GAL_5 -sn GAL_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group12_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group13 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GEO_2 -sn GEO_6 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group13_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group14 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn GRE_1 -sn GRE_2 -sn GRE_3 -sn GRE_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group14_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group15 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group15_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group16 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LNL_3 -sn MEL_1 -sn SKI_5 -sn ROT_13 -sn SPU_8 -sn ERS_1 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group16_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group17 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn LOS_1 -sn LOS_2 -sn LOS_6 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group17_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group18 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group18_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group19 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn OSW_2 -sn OSW_3B -sn OSW_6 -sn OSW_7 -sn OSW_8 -sn OSW_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group19_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group20 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn RUN_1 -sn RUN_4 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group20_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group21 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn SCU_15 -sn SCU_16 -sn SCU_19 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group21_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group22 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn SKN_1 -sn SKN_2 -sn SKN_5 -sn SKN_8 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group22_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group23 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn TNG_3 -sn TNG_5 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group23_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group24 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn TRO_1 -sn TRO_5 -sn TRO_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group24_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group25 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn UIG_1 -sn UIG_10 -sn UIG_2 -sn UIG_3 -sn UIG_4 -sn UIG_6 -sn UIG_7 -sn UIG_8 -sn UIG_9 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group25_ALT.fa
rm imptop-rated-overall_ALT.vcf*

#Group26 - filter vcf and generate ALT sequence
gatk SelectVariants --restrict-alleles-to BIALLELIC -V $vcf --select "AF > 0.5" -sn VAG_4 -sn VAG_7 -L Cexcelsa_scaf_2:4113454-4117809 -O imptop-rated-overall_ALT.vcf.gz
samtools faidx $ref Cexcelsa_scaf_2:4113454-4114535 Cexcelsa_scaf_2:4117100-4117330 Cexcelsa_scaf_2:4117593-4117809 | bcftools consensus imptop-rated-overall_ALT.vcf.gz -o imp_top_rated_SNP/top-rated-group26_ALT.fa
rm imptop-rated-overall_ALT.vcf*
