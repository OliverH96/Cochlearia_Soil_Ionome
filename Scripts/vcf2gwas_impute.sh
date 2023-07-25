#!/bin/bash
#SBATCH --job-name=vcf2gwas_imp
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64g
#SBATCH --time=150:00:00
#SBATCH --output=/gpfs01/home/mbxoh2/OandE/%x.out
#SBATCH --error=/gpfs01/home/mbxoh2/OandE/%x.err

source $HOME/.bash_profile
conda activate vcf2gwas

cd /gpfs01/home/mbxoh2/genome_data/

#set variables
lists=./cochlearia/Element_lists/*imputed.csv
gff=./reference/C_excelsa_V5_braker2_wRseq.gff3
vcf=./all_bi_best_snps.vcf.gz
out=./cochlearia/gwas_out/

# Redirect output to a log file
log_file="/gpfs01/home/mbxoh2/OandE/vcf2gwas_imp_log.txt"
exec &> >(tee -a "$log_file")

#run Gwas for each element phenotype
for list in $lists
do
        name=$(basename "$list" .csv)
        mkdir "$out""$name"
        vcf2gwas -v $vcf -pf $list -gf $gff -ap -lm 4 -o "$out""$name"
done

