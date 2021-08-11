# Rank signif exclude two

## Ranking of all 26 rare disease groups, update for revision, includes
#' 1. co-expression: core
#' 2. different PPIs (large/small scales)
#' The script is to compare the performances
## Ize Buphamalai
## June 2021
## This uses the LCC computation from all genes -> specifically exclude each layers to compare performances


library(igraph)
library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
distinct <- dplyr::distinct

source("../functions/process_rank_network.R")

# 1. Load disease network prior matrix ------------
## load the LCC results and prepare data for the original prior matrix
result_folder = "../cache/output/Orphageneset_rare/"
result_df = readRDS(paste0(result_folder, "LCC_and_distance_calculation_results.RDS"))


## load the LCC results and prepare data for the newly added prior matrix
result_df_additional_net = readRDS("../cache/output/Orphageneset_rare_additionalNetworks/LCC_and_distance_calculation_results.RDS")

result_df_combn <- rbind(result_df, 
                         result_df_additional_net %>% filter(network %in% c("coex_core")))

## from LCC_analysis_new: process the LCC
result_df_combn = process_LCC_result(result_df_combn, network_annotate = F)

### extract the most significant networks for all diseases
all_LCC_val_signif = result_df_combn %>% dplyr::filter(correctedPval < 0.05) %>% dplyr::select(name, network, LCC.zscore)



# 2. Load disease gene list ------------
Orphanet_df = process_disease_genes_data("../data/table_disease_gene_assoc_orphanet_genetic.tsv", min_gene = 20, max_gene = 2000)
Orphanet_df = Orphanet_df$disgene_list


# 3. Load edgelists ------------
network_dirs = "../data/network_edgelist_combn/"

#all_networks = all_LCC_val_signif %>% pull(network) %>% unique
all_networks =  str_remove(list.files(network_dirs), pattern = ".tsv")


el_all = list()
for(i in all_networks){
  el_all[[i]] = read_tsv(paste0(network_dirs, i,".tsv"), col_names = c("A","B"), skip = 1, col_types = 'cc')
}

el_all = lapply(el_all, process_edgelist)


## train and test sets
all_diseases = names(Orphanet_df)

# 4 Perform rank retrieval for different network sets ------------
rank_networks <- list()


rank_networks[["Signif_withoutPPI"]] <-process_rank_network(all_diseases = all_diseases,
                                                             el_all = el_all, 
                                                             weighted = T, 
                                                             disease_specific = T, 
                                                             to_exclude = c("ppi"))

rank_networks[["Signif_withoutGOBP"]] <-process_rank_network(all_diseases = all_diseases,
                                                           el_all = el_all, 
                                                           weighted = T, 
                                                           disease_specific = T, 
                                                           to_exclude = c("GOBP"))


rank_networks[["Signif_withoutGOMF"]] <-process_rank_network(all_diseases = all_diseases,
                                                           el_all = el_all, 
                                                           weighted = T, 
                                                           disease_specific = T, 
                                                           to_exclude = c("GOMF"))


rank_networks[["Signif_withoutHP"]] <-process_rank_network(all_diseases = all_diseases,
                                                           el_all = el_all, 
                                                           weighted = T, 
                                                           disease_specific = T, 
                                                           to_exclude = c("HP"))


rank_networks[["Signif_withoutMP"]] <-process_rank_network(all_diseases = all_diseases,
                                                           el_all = el_all, 
                                                           weighted = T, 
                                                           disease_specific = T, 
                                                           to_exclude = c("MP"))


rank_networks[["Signif_withoutCoPathway"]] <-process_rank_network(all_diseases = all_diseases,
                                                           el_all = el_all, 
                                                           weighted = T, 
                                                           disease_specific = T, 
                                                           to_exclude = c("reactome_copathway"))

#el_all_mod = el_all
#el_all_mod$ppi =  process_edgelist(read_tsv("../data/network_edgelist_additional/PPI_HIPPIELargeScale.tsv", col_names = c("A","B"), skip = 1, col_types = 'cc'))
#rank_networks[["Signif_withoutLargeScale"]] <-process_rank_network(all_diseases = all_diseases,
#                                                                  el_all = el_all, 
#                                                                  weighted = T, 
#                                                                  disease_specific = T, 
#                                                                  to_exclude = c( "MP", "HP", "GOMF","GOBP", "reactome_copathway"))


# Process fold results ###########################

rank_networks_processed <- lapply(rank_networks, function(x) 
  lapply(x, function(df) 
    lapply(1:10, function(iter) 
      df %>% filter(fold == iter) %>% select(trueset, rank))))


## compute AUC from the rank through CVAUC - combining cross validation results
pacman::p_load(cvAUC)
aucvals <- pbapply::pblapply(rank_networks_processed,  function(x) lapply(x, auc_output))


######
AUC_folds_each <- lapply(1:length(aucvals), function(x) auc_to_df(aucvals[[x]], label = names(aucvals)[x]))

AUC_folds <- bind_rows(AUC_folds_each)

saveRDS(AUC_folds, "../cache/fold_cv_processed_results_revision_excludesSingleLayers.RDS")


#############
# Merge GOs and HP & MPs

# rank_networks <- list()
# 
# rank_networks[["Signif_withoutPhenotype"]] <-process_rank_network(all_diseases = all_diseases,
#                                                                   el_all = el_all, 
#                                                                   weighted = T, 
#                                                                   disease_specific = T, 
#                                                                   to_exclude = c("HP", "MP"))
# 
# rank_networks[["Signif_withoutGO"]] <-process_rank_network(all_diseases = all_diseases,
#                                                            el_all = el_all, 
#                                                            weighted = T, 
#                                                            disease_specific = T, 
#                                                            to_exclude = c("GOMF", "GOBP"))


# Process fold results ###########################

# rank_networks_processed <- lapply(rank_networks, function(x) 
#   lapply(x, function(df) 
#     lapply(1:10, function(iter) 
#       df %>% filter(fold == iter) %>% select(trueset, rank))))
# 
# 
# ## compute AUC from the rank through CVAUC - combining cross validation results
# pacman::p_load(cvAUC)
# aucvals <- pbapply::pblapply(rank_networks_processed,  function(x) lapply(x, auc_output))
# 
# 
# ######
# AUC_folds_each <- lapply(1:length(aucvals), function(x) auc_to_df(aucvals[[x]], label = names(aucvals)[x]))
# 
# AUC_folds <- bind_rows(AUC_folds_each)
# 
# saveRDS(AUC_folds, "../cache/fold_cv_processed_results_revision_excludesOntology_merged.RDS")