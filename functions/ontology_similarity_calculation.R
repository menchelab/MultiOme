# Compute the semantic similarity among disease genes


# set up and load packages
##############
source("../functions/ontology_similarity_functions.R")
source("../functions/readdata_functions.R")

## read disease gene-disease relationship file
##############
disgene = process_disease_genes_data("../data/table_disease_gene_assoc_orphanet.tsv")
disgene_pass_df = disgene$disgene_df
disgene_list = disgene$disgene_list

# restrict the disease list that pass the criteria
disgene_pass_list = disgene_list[disgene_pass_df$name]

# all the table in the database to use
table_list = c("gobp", "gomf", "mp", "hp", "go_minfreq")

# main calculation - take time
term_similarity_df = diseaseset_sim_measures(table_list, disgene_pass_list)

# save results
write_tsv(term_similarity_df, path = "./output/disease_gene_ontology_similarity.tsv")
