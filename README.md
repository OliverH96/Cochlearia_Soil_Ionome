# Cochlearia_Soil_Ionome
This repository contains data and code produced and used in the MSc Bioinformatics final individual research project of Oliver Hargreaves, supervised by Sian Bray. 

## Soil Ionome Imputation
For soil ionome data imputation, gridding and imputation was performed using [Surfer®](https://www.goldensoftware.com/products/surfer)

Grids were generated using manually-fitted variograms (Gridding/MFV/) and auto-fitted variograms (Gridding/AFV/). Grids for each element can be found in their respective directories. Gridding/MFV/imputation_setup.xlsx contains the variables used to fit MFV grids.


Validation correlation plots were generated using the correlation_plots.py and correlation_plots_uneditedVG.py scripts. FOREGS imputed data used for validation is contained within Gridding/FOREGS_Impute_raw.xlsx

## Genomic Information Field Theory (GIFT)
Bash scripts used call gift.py to perfrom GIFT analysis with non-imputed and imputed datasets can be found in Scripts/GIFT_no_impute and Scripts/GIFT_impute respectively. Note, GIFT original code is currently unpublished and was used with permission of its authors.

Outlier SNPs were annotated using Extract_function.py to match them with any gene regions in which they occured as well as to provide a brief functional description of any known *Arabidopsis thaliana* orthologs.
