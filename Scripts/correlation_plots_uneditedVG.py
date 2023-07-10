#Import modules
import os

#Elements list and change directory
elements=["NA2O", "MGO", "P2O5", "CAO", "MNO", "RB", "SR", "FE2O3", "K2O", "ZN"]
act="_T_Act"
imp=" No Edit Variogram Imp"
os.chdir("C:/Users/Oliver/Desktop/UoN/Individual_project/Imputation_Grids/Validation")

#For loop to produce correlation plots for imputed and actual data
for i in elements:
     if i == "NA2O":
          name='"Na"[2]*"O"~"'
     elif i == "MGO":
          name='"MgO '
     elif i == "P2O5":
          name='"P"[2]*"O"[5]~"'
     elif i == "K2O":
          name='"K"[2]*"O"~"'
     elif i == "CAO":
          name='"CaO '
     elif i == "MNO":
          name='"MnO '
     elif i == "FE2O3":
          name='"Fe"[2]*"O"[3]~"'
     elif i == "ZN":
          name='"Zn '
     elif i == "RB":
          name='"Rb '
     elif i == "SR":
          name='"Sr '
     out_file3=open('plot.r', 'w')
     out_file3.write(f'library(ggpubr)\n'+
          f'library(readxl)\n'+
          f'foregsdata <- read_excel("C:/Users/Oliver/Desktop/UoN/Individual_project/Imputation_Grids/Validation/FOREGS_Impute_raw.xlsx")\n'+
          f'ggscatter(foregsdata, x="{i}{imp}", y="{i}{act}", add = "reg.line", conf.int = TRUE) +\n'+
          f'  ylab({name}FOREGS Actual Data") + xlab({name}Imputed Data") +\n'+
          f'  stat_cor(method = "spearman")\n'+
          f'ggsave("{i}NoEditVGplot.tiff", device = "tiff")\n')
     out_file3.close()
     
     os.system('Rscript plot.r')
