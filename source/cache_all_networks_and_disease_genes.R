# cache all network layers for faster loading, and disease association

# location where cache is installed
loc <- "../cache/all_networks_and_disease_genes_cached.RData"


if(!file.exists(loc)){
  source("../functions/readdata_functions.R")
  source("../functions/readdata_functions.R")
  source("../functions/process_LCC_result.R")
  
  pacman::p_load(treemap, igraph)
  
  g = process_graph_data("../data/network_edgelist_combn//")
  
  rare_genetic_result_folder = "../cache/output/Orphageneset_rare/"
  ## process the results LCC significance
  result_folder = "../cache/output/Orphageneset_rare/"
  result_df = readRDS(paste0(result_folder, "LCC_and_distance_calculation_results.RDS"))
  
  
  ## load the LCC results and prepare data for the newly added prior matrix
  result_df_additional_net = readRDS("../cache/output/Orphageneset_rare_additionalNetworks/LCC_and_distance_calculation_results.RDS")
  
  result_df <- rbind(result_df, 
                           result_df_additional_net %>% filter(network %in% c("coex_core")))
  
  
  rare_genetic_diseases_genes <- process_disease_genes_data("../data/table_disease_gene_assoc_orphanet_genetic.tsv", 20,2000)
 
  ## from LCC_analysis_new: process the LCC
  processed_result_df = process_LCC_result(result_df, network_annotate = T) %>%
    filter(name %in% names(rare_genetic_diseases_genes$disgene_list))
  
   
  individual_diseases_genes <- read_tsv("../data/orphaNet_disease_gene_association_with_roots.tsv")
  
  save(list = c("g", "processed_result_df", "rare_genetic_diseases_genes", "individual_diseases_genes"), 
       file = loc)
} else{
  load(loc)
}