# load topGO
library("topGO")

setwd("C://Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/no_impute_out/Na/")
# Open my files

# This function converts the file to the gene2GO format that topGO requires
go <- readMappings(file="C://Users/Oliver/Desktop/UoN/Individual_project/Cochlearia_Thaliana_GO_universe_restrictive_id2gos_no_obsolete_nonredundant.tsv", 
                   sep = "\t", IDsep = ",")

# All the genes in my scan (i.e. all the genes in regions from the bpm file)
genesInScan <- read.csv("C://Users/Oliver/Desktop/UoN/Individual_project/Cochlearia_all_genes.txt", 
                        header = FALSE)

# My genes of interest (i.e. Fst top 1%)
top1perc <- read.csv("C://Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/no_impute_out/Na/genesofinterest.txt", 
                     header = FALSE)

# Create a factor where the indecies are the gene names and the factors are 1/0 for "candidate gene"/"scan gene but not candidate gene"

# First convert top1perc from a dataframe to a vector
top1perc_vector <- as.vector(top1perc[,'V1'])

# Then convert genesInScan to a vector
genesInScan_vector <- as.vector(genesInScan[,'V1'])

# Then make the factor
allGenes_factor <- factor(as.integer(genesInScan_vector %in% top1perc_vector))
names(allGenes_factor) = genesInScan_vector

# Website used for help: https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
# Make the topGOdata object
go_data_BP <- new("topGOdata", 
               ontology="BP", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene

go_data_CC <- new("topGOdata", 
               ontology="CC", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene

go_data_MF <- new("topGOdata", 
               ontology="MF", 
               allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan
               annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"
               gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms
               nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene


# Run the GO analysis classic fisher
res_BP_classicFisher <- runTest(go_data_BP, statistic = "fisher")
res_CC_classicFisher <- runTest(go_data_CC, statistic = "fisher")
res_MF_classicFisher <- runTest(go_data_MF, statistic = "fisher")

# Run the GO analysis conserative elim fisher
resBP_elimFisher <- runTest(go_data_BP, algorithm = "elim", statistic = "fisher")
resCC_elimFisher <- runTest(go_data_CC, algorithm = "elim", statistic = "fisher")
resMF_elimFisher <- runTest(go_data_MF, algorithm = "elim", statistic = "fisher")

# Check it
res_BP_classicFisher
res_CC_classicFisher
res_MF_classicFisher

# Check it
resBP_elimFisher
resCC_elimFisher
resMF_elimFisher

# Show combined table
allGO_BP=usedGO(go_data_BP)
allGO_CC=usedGO(go_data_CC)
allGO_MF=usedGO(go_data_MF)

allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,
                     elimFisher = resBP_elimFisher,
                     orderBy = "elimFisher",
                     topNodes=length(allGO_BP),
                     numChar = 250)

allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,
                     elimFisher = resCC_elimFisher,
                     orderBy = "elimFisher",
                     topNodes=length(allGO_CC),
                     numChar = 250)

allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,
                     elimFisher = resMF_elimFisher,
                     orderBy = "elimFisher",
                     topNodes=length(allGO_MF),
                     numChar = 250)

# Write the combined tables to files
write.csv(allResBP, "allResBP.csv")
write.csv(allResCC, "allResCC.csv")
write.csv(allResMF, "allResMF.csv")




