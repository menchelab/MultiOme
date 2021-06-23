## Ranking of all 26 rare disease groups, update for revision, includes
#' 1. co-expression: core
#' 2. different PPIs (large/small scales)
#' The script is to compare the performances
## Ize Buphamalai
## June 2021
## This uses the LCC computation from all genes


library(igraph)
source("../functions/process_rank_network.R")

# 1. Load disease network prior matrix ------------
## load the LCC results and prepare data for the original prior matrix
result_folder = "../cache/output/Orphageneset_rare_CoExwithCore/"
result_df = readRDS(paste0(result_folder, "LCC_and_distance_calculation_results.RDS"))


## load the LCC results and prepare data for the newly added prior matrix
result_df_additional_net = readRDS("../cache/output/Orphageneset_rare_additionalNetworks/LCC_and_distance_calculation_results.RDS")

result_df_combn <- rbind(result_df, 
                   result_df_additional_net %>% filter(network %in% c("PPI_HIPPIECurated" , "PPI_HIPPIELargeScale" ,"coex_core")))

## from LCC_analysis_new: process the LCC
result_df_combn = process_LCC_result(result_df_combn)

### extract the most significant networks for all diseases
all_LCC_val_signif = result_df_combn %>% dplyr::filter(correctedPval < 0.05) %>% dplyr::select(name, network, LCC.zscore)



# 2. Load disease gene list ------------
Orphanet_df = process_disease_genes_data("../data/table_disease_gene_assoc_orphanet_genetic.tsv", min_gene = 20, max_gene = 2000)
Orphanet_df = Orphanet_df$disgene_list


# 3. Load edgelists ------------
network_dirs = "../data/network_edgelist_combn/"

#all_networks = all_LCC_val_signif %>% pull(network) %>% unique
all_networks = str_remove(list.files(network_dirs), pattern = ".tsv")


el_all = list()
for(i in all_networks){
  el_all[[i]] = read_tsv(paste0(network_dirs, i,".tsv"), col_names = c("A","B"), skip = 1, col_types = 'cc')
}

el_all = lapply(el_all, process_edgelist)


## train and test sets
all_diseases = names(Orphanet_df)

# 4 Perform rank retrieval for different network sets ------------

rank_networks <- list()
## 1) Ranking for PPI Large scale only

rank_networks[["PPI_largescale"]] <-process_rank_network(network_set = "PPI_HIPPIELargeScale", 
                                                         all_diseases = all_diseases,
                                                         el_all = el_all, 
                                                         weighted = F, disease_specific = F, to_exclude = NULL)

## 2) Ranking for PPI curated only
rank_networks[["PPI_curated"]] <-process_rank_network(network_set = "PPI_HIPPIECurated", 
                                                         all_diseases = all_diseases,
                                                         el_all = el_all, 
                                                         weighted = F, disease_specific = F, to_exclude = NULL)

## 3) signif - upper - ppi + ppi_large 
upper_layers <- c("GOBP","GOMF","reactome_copathway","HP","MP")

rank_networks[["Largescale_signif"]] <- process_rank_network(all_diseases = all_diseases,
                                                               el_all = el_all, 
                                                               weighted = F, disease_specific = T, 
                                                               to_exclude = c(upper_layers, "ppi", "PPI_HIPPIECurated"))
  

## 4)  all layers (updated by adding core coexpression)

all_network_df <- tibble(network = setdiff(unique(all_LCC_val_signif$network), c("PPI_HIPPIELargeScale", "PPI_HIPPIECurated")), LCC.zscore = 1)


rank_networks[["all"]] <-  process_rank_network(network_weight_df_custom=all_network_df,
                                                all_diseases = all_diseases,
                                                el_all = el_all, 
                                                weighted = F, 
                                                disease_specific = F, 
                                                to_exclude = NULL)


## 5)  Signif (include coex-core): update only those where coex_core is significant
diseases_to_consider <- all_LCC_val_signif %>% filter(grepl("core", network)) %>% pull(name)


rank_networks[["signif_withCoexCore"]] <- process_rank_network(all_diseases = diseases_to_consider,
                                                             el_all = el_all, 
                                                             weighted = T, 
                                                             disease_specific = T, 
                                                             to_exclude = c("PPI_HIPPIELargeScale", "PPI_HIPPIECurated"))


## 6) Signif - upper layers (like #3) but weighted
upper_layers <- c("GOBP","GOMF","reactome_copathway","HP","MP")

rank_networks[["Largescale_signif_weighted"]] <- process_rank_network(all_diseases = all_diseases,
                                                             el_all = el_all, 
                                                             weighted = T, 
                                                             disease_specific = T, 
                                                             to_exclude = c(upper_layers, "ppi", "PPI_HIPPIECurated"))


## 7) scale_specific
# 7.1 - co-essential alone
rank_networks[["co-essential"]] <- process_rank_network(network_set = "co-essential", 
                                                        all_diseases = all_diseases,
                                                        el_all = el_all, 
                                                        weighted = F, disease_specific = F, to_exclude = NULL)

# 7.2 - reactome alone
rank_networks[["reactome_copathway"]] <- process_rank_network(network_set = "reactome_copathway", 
                                                        all_diseases = all_diseases,
                                                        el_all = el_all, 
                                                        weighted = F, disease_specific = F, to_exclude = NULL)

# 7.3 - phenotype alone
phenotype_weight_df <- tibble(network = c("HP","MP"), LCC.zscore = 1)

rank_networks[["phenotypes"]] <- process_rank_network(network_weight_df_custom=phenotype_weight_df,
                                                      all_diseases = all_diseases,
                                                      el_all = el_all, 
                                                      weighted = F, 
                                                      disease_specific = F, 
                                                      to_exclude = NULL)

# 7.4 - co-expression alone - only consider the significant ones (tissues)
rank_networks[["coex_Signif"]] <- process_rank_network(all_diseases = all_diseases,
                                                        el_all = el_all, 
                                                        weighted = T, 
                                                        disease_specific = T, 
                                                        to_exclude = names(el_all)[!grepl("coex", names(el_all))])

# 7.5 - GO networks
go_weight_df <- tibble(network = c("GOBP","GOMF"), LCC.zscore = 1)

rank_networks[["GO"]] <-  process_rank_network(network_weight_df_custom=go_weight_df,
                                               all_diseases = all_diseases,
                                               el_all = el_all, 
                                               weighted = F, 
                                               disease_specific = F, 
                                               to_exclude = NULL)

# 7.6 - core co-expression
rank_networks[["coex_core"]] <- process_rank_network(network_set = "coex_core", 
                                                              all_diseases = all_diseases,
                                                              el_all = el_all, 
                                                              weighted = F,disease_specific = F, to_exclude = NULL)

# Process fold results ###########################

rank_networks_processed <- lapply(rank_networks, function(x) 
  lapply(x, function(df) 
    lapply(1:10, function(iter) 
      df %>% filter(fold == iter) %>% select(trueset, rank))))

saveRDS(rank_networks_processed, "../cache/rank_10foldCV_revised_26_diseases.RDS")


## compute AUC from the rank through CVAUC - combining cross validation results
aucvals <- pbapply::pblapply(rank_networks_processed,  function(x) lapply(x, auc_output))


######
AUC_folds_each <- lapply(1:length(aucvals), function(x) auc_to_df(aucvals[[x]], label = names(aucvals)[x]))

AUC_folds <- bind_rows(AUC_folds_each)

saveRDS(AUC_folds, "../cache/fold_cv_processed_results_revision.RDS")



