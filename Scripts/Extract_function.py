#A python script to  obtain gene IDs and functional descriptions for top SNPs from GIFT using a gff and function file
#Example command:
#python /c/Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/Extract_function.py -c /c/Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/impute_out/Na/top_output.csv -g /c/Users/Oliver/Desktop/UoN/Individual_project/Reference/C_excelsa_V5_braker2_wRseq.gff3 -f /c/Users/Oliver/Desktop/UoN/Individual_project/Reference/1-2-1_hits_all_gene_descriptions.tsv -o /c/Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/impute_out/Na

import pandas as pd
import argparse, os

#Set parser arguments for directory input
parser = argparse.ArgumentParser()
parser.add_argument('-c', type=str, metavar='csv_path', required=True, help='Path to input csv containing top SNPs')
parser.add_argument('-g', type=str, metavar='gff_path', required=True, help='Path to gff file')
parser.add_argument('-f', type=str, metavar='function_path', required=True, help='Path to gene function file')
parser.add_argument('-o', type=str, metavar='dir_path', required=True, help='Output directory path to save annotated gene list')
args = parser.parse_args()

#read SNP and gff files as dataframes
gff=pd.read_csv(args.g, sep='\t', header=None)
snp=pd.read_csv(args.c)
function=pd.read_csv(args.f, sep='\t', header=None)
os.chdir(args.o)
print(function)

#Keep only gene regions in gff
filtered_rows=[]

for index, row in gff.iterrows():
    if 'gene' in row.values:
        filtered_rows.append(row.to_dict())
filtered_gff=pd.DataFrame(filtered_rows)

#Obtain start and end regions of each gene as well as the scaffold
final_filtered_rows=[]
for index, row in filtered_gff.iterrows():
    chr_value=row[0]
    start=row[3]
    end=row[4]
    
    #Match snp Chromosome values with genes from gff on the same chromosome
    snp_filtered=snp[snp['CHR'] == chr_value]
    
    #Check these potential matches to see whether snp BP falls within gene region
    for snp_index, snp_row in snp_filtered.iterrows():
        bp_value = snp_row['BP']
        if start <= bp_value <= end:
            final_filtered_rows.append([chr_value, bp_value, row[8], row[2]])
            
            
#Create intermediate dataframes
df=pd.DataFrame(final_filtered_rows, columns=['CHR', 'BP', 'ID', 'Type'])
df['ID'] = df['ID'].str.replace('ID=', '').str.replace(';', '.')

#Create new column to describe gene function and A. Thaliana Ortholog ID
df['Ortholog ID'] =''
df['Function']=''

#Find partial matches between ID in intermediate dataframe and gene IDs in gene function file:
for index, row in df.iterrows():
    id_value=row['ID']
    #Check for partial matches in function file
    matches = function[function.iloc[:, 0].str.contains(id_value, na=False)]
    #If there is a match, extract the functional description
    if not matches.empty:
        df.at[index, 'Ortholog ID'] = matches.iloc[0, 1]
        df.at[index, 'Function'] = matches.iloc[0, 3]


#Save final dataframe to csv       
print(df.shape[0],'SNPs found within gene regions. Saving to csv.')
df.to_csv('Annotated-SNPs.csv', index=False)
