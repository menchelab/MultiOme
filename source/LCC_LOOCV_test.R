#### Multiscale RWR for all the Orphanet: using leave oue out cross validation -> to simulate scenarios for unknown genes
## Ize Buphamalai
## 15 Feb 2021
## This uses the LCC computation from all genes

source("../functions/readdata_functions.R")
source("../functions/LCC_functions.R")
source("../functions/process_edgelist.R")
source("../functions/process_LCC_result.R")
source("../functions/weighted_multiplex_propagation.R")
source("../functions/RWR_get_allnodes.R")
source("../functions/RWR_transitional_matrix.R")
source("../functions//RWR.R")
library(igraph)


## load the LCC results 

result_folder = "../cache/output/Orphageneset_rare/"
result_df = readRDS(paste0(result_folder, "LCC_and_distance_calculation_results.RDS"))

## from LCC_analysis_new: process the LCC
result_df = process_LCC_result(result_df)

### extract the most significant networks for all diseases
all_LCC_val_signif = result_df %>% dplyr::filter(correctedPval < 0.05) %>% dplyr::select(name, network, LCC.zscore)

## load the gene sets
Orphanet_df = process_disease_genes_data("../data/table_disease_gene_assoc_orphanet_genetic.tsv", min_gene = 20, max_gene = 2000)
Orphanet_df = Orphanet_df$disgene_list


## read all edgelists for all significant networks
network_dirs = "../data/network_edgelists/"

all_networks = all_LCC_val_signif %>% pull(network) %>% unique

el_all = list()
for(i in all_networks){
  el_all[[i]] = read_tsv(paste0(network_dirs, i,".tsv"), col_names = c("A","B"), skip = 1)
}

el_all = lapply(el_all, process_edgelist)


## train and test sets
all_diseases = names(Orphanet_df)



## Ranking for significant networks

## to store ranking results
rank_df_all_folds = list()

for(disease_current in all_diseases){
  
  print(disease_current)
  
  # set seed
  disease_genes = Orphanet_df[[disease_current]]
  n = length(disease_genes)
  
  
  LCC_val_signif = all_LCC_val_signif %>% dplyr::filter(name == disease_current) %>% select(network, LCC.zscore)
  print(sprintf("Number of network taken into account is %i", nrow(LCC_val_signif)))
  
  el = el_all[LCC_val_signif$network]
  gene_allnet = get_allnodes(el)
  supraadj = supraadjacency_compute(LCC_val_signif, el = el)

  
  # each retrievat sets controlled by parameter i, i = 1:k
  rank_df = list()
  
  for(i in 1:n){
    seed_sets = disease_genes[-i]
    retrieval_sets = disease_genes[i]
    
    ## perform retrieval
    
    ranks <- weighted_multiplex_propagation(
      seedset = seed_sets, 
      trueset = disease_genes, 
      weighted_layer_df = LCC_val_signif,
      network_dirs = network_dirs, 
      el = el, 
      gene_allnet = gene_allnet, 
      Stest = supraadj)
 
    
    rank_df[[i]] = ranks %>% dplyr::filter(trueset)
    
     }
  
 
  all_ranks <-  bind_rows(rank_df, .id = "fold") %>% mutate(max = length(gene_allnet))
  rank_df_all_folds[[disease_current]] = all_ranks
  
}

#})

saveRDS(rank_df_all_folds, "../cache/LCC_LOOCV_ranking_rare_genetic_diseases.RDS")

##############

# prepare edge list, and supra-adjacency matrix (PPI alone)

LCC_val_signif = tibble(network = "ppi", LCC.zscore = 1) # just to add PPI

el = el_all[LCC_val_signif$network]
gene_allnet = get_allnodes(el)
supraadj = supraadjacency_compute(LCC_val_signif, el = el)


rank_df_all_folds_only_PPI = list()
for(disease_current in all_diseases){
  
  print(disease_current)
  
  # set seed
  disease_genes = Orphanet_df[[disease_current]]
  n = length(disease_genes)
  
  
  rank_df = list()
  for(i in 1:n){
    seed_sets = disease_genes[-i]
    retrieval_sets = disease_genes[i]   
    ## perform retrieval
    
    ranks = weighted_multiplex_propagation(
      seedset = seed_sets, 
      trueset = disease_genes, 
      weighted_layer_df = LCC_val_signif,
      network_dirs = network_dirs, 
      el = el, 
      gene_allnet = gene_allnet, 
      Stest = supraadj)
    
    rank_df[[i]] = ranks %>% dplyr::filter(trueset)
  }
  
  all_ranks<-  bind_rows(rank_df, .id = "fold") %>% mutate(max = length(gene_allnet))
  rank_df_all_folds_only_PPI[[disease_current]] = all_ranks
  
}





saveRDS(rank_df_all_folds_only_PPI, "../cache/LCC_LOOCV_ranking_rare_genetic_diseases_only_PPI.RDS")


##############

# prepare edge list, and supra-adjacency matrix (all networks)

LCC_val_signif = tibble(network = unique(result_df$network), LCC.zscore = 1) # all networks

el = el_all[LCC_val_signif$network]
gene_allnet = get_allnodes(el)
supraadj = supraadjacency_compute(LCC_val_signif, el = el)


rank_df_all_folds_allnets = list()
for(disease_current in all_diseases){
  
  print(disease_current)
  
  # set seed
  disease_genes = Orphanet_df[[disease_current]]
  n = length(disease_genes)
  
  
  rank_df = list()
  for(i in 1:n){
    seed_sets = disease_genes[-i]
    retrieval_sets = disease_genes[i]   
    ## perform retrieval
    
    ranks = weighted_multiplex_propagation(
      seedset = seed_sets, 
      trueset = disease_genes, 
      weighted_layer_df = LCC_val_signif,
      network_dirs = network_dirs, 
      el = el, 
      gene_allnet = gene_allnet, 
      Stest = supraadj)
    
    rank_df[[i]] = ranks %>% dplyr::filter(trueset)
  }
  
  all_ranks<-  bind_rows(rank_df, .id = "fold") %>% mutate(max = length(gene_allnet))
  rank_df_all_folds_allnets[[disease_current]] = all_ranks
  
}


saveRDS(rank_df_all_folds_allnets, "../cache/LCC_LOOCV_ranking_rare_genetic_diseases_allnets.RDS")
