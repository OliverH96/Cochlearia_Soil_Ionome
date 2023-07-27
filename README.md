# Cochlearia_Soil_Ionome
This repository contains data and code produced and used in the MSc Bioinformatics final individual research project of Oliver Hargreaves, supervised by Sian Bray. 

## Soil Ionome Imputation
For soil ionome data imputation, gridding and imputation was performed using [Surfer®](https://www.goldensoftware.com/products/surfer)

Grids were generated using manually-fitted variograms (Gridding/MFV/) and auto-fitted variograms (Gridding/AFV/). Grids for each element can be found in their respective directories. Gridding/MFV/imputation_setup.xlsx contains the variables used to fit MFV grids.


Validation correlation plots were generated using the correlation_plots.py and correlation_plots_uneditedVG.py scripts. FOREGS imputed data used for validation is contained within Gridding/Validation/FOREGS_Impute_raw.xlsx as well as correlation plots for the 10 element validation set MFV and AFV.

Actual and imputed ionome data can be found in the Cochlearia_Soil_Master csv files.

## GWAS and Genomic Information Field Theory (GIFT)
Despite not vcf2gwas not finding any significant genotype/phenotype associations, bash scripts used to run vcf2gwas can be found in Scripts/

Bash scripts used call gift.py to perfrom [GIFT](https://iopscience.iop.org/article/10.1088/1478-3975/ac99b3) analysis with non-imputed and imputed datasets can be found in Scripts/GIFT_no_impute and Scripts/GIFT_impute respectively. Note, GIFT original code is currently unpublished and was used with permission of its authors.

Outlier SNPs were annotated using Extract_function.py to match them with any gene regions in which they occured as well as to provide a brief functional description of any known *Arabidopsis thaliana* orthologs.

GIFT output manhattan plots for each element are located in GIFT_Output directory. Full output csv are not present here due to file size. The analysed top 0.01% SNPs for Na from both datasets are included, along with scripts and full output for GO and STRING analyses. Cochlearia GO term mappings an be found in Cochlearia_Thaliana_GO_universe_restrictive_id2gos_no_obsolete_nonredundant.tsv

## Consensus Sequences and Predicted Structures
Scripts used to generate group consensus sequences can be found in Scripts/Consensus_Sequences
Groupings and vcf files used for consensus sequence generation can be found in Consensus_Sequences
Final, EMBOSS generated g40302.t1 reference and alternate consensus sequences can be found in Consensus_Sequences/ alongside translated sequences and predicted structural models. All sequences are in fasta format, models in pdb format.

Multiple sequence alignment of HKT1 orthologs used in observations of conserved residues and motifs can be found in HKT1_Ortholog_Alignment.clustal_num
