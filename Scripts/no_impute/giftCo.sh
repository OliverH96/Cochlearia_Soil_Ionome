#!/bin/bash
#SBATCH --job-name=giftCo
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=179g
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/mbxoh2/OandE/%x.out
#SBATCH --error=/gpfs01/home/mbxoh2/OandE/%x.err

source $HOME/.bash_profile
conda activate gift

cd /gpfs01/home/mbxoh2/genome_data/cochlearia/gift/no_impute_out

# Set variables
list=/gpfs01/home/mbxoh2/genome_data/cochlearia/Element_lists/Cochlearia_Soil_Master.csv
vcf=/gpfs01/home/mbxoh2/genome_data/cochlearia/gift/all_bi_best_snps.vcf
mkdir Co
out=/gpfs01/home/mbxoh2/genome_data/cochlearia/gift/no_impute_out/Co/
script=/gpfs01/home/mbxoh2/genome_data/cochlearia/gift/gift.py


# Run Gwas for each element phenotype
python3 -u $script -v "$vcf" -f "$list" -p Co -o "$out"output.csv >> "$out"log.txt
