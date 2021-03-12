# The file to read and process immediate orphanet gene-disease association
# Ize Buphamalai
# Last update: 09.04.2020

############
# read orphanet data
############

orphaGeneDisease = read_csv("../Disease_gene_assoc/Orphanet/output/processsed_Orphanet_genedisease_association.csv", col_names = c("disease", "genes"), col_types = 'dc')
orphaGeneDisease$genes = sapply(orphaGeneDisease$genes, function(x) strsplit(x, ";")[[1]])
orphaGeneDisease$n_genes = sapply(orphaGeneDisease$genes, length)
orphaGeneDisease$orphaID = paste0("Orphanet:", orphaGeneDisease$disease)

# convert to a list view
orphaGeneDiseaseList = orphaGeneDisease$genes
names(orphaGeneDiseaseList) = orphaGeneDisease$orphaID





