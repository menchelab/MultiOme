# sample the PPI to the same amount of rows as the large scale PPI

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
el_PPI = process_edgelist(read_tsv("../data/network_edgelist_additional/ppi_HuRI.tsv", col_names = c("A","B"), skip = 1, col_types = 'cc'))


## train and test sets
all_diseases = names(Orphanet_df)

# 4 Perform rank retrieval for different network sets ------------
# randomly take 
evals <- 50304
#split(d, ceiling(seq_along(d)/20))

rank_networks_PPI <- list()

for(i in 1:10){
  print(i)
  set.seed(i)
#  label <- paste(i)
  to_select <- sample(1:nrow(el_PPI), evals, replace = F)
  el_used <- el_PPI[to_select,]
  g_used <- graph_from_data_frame(el_used)
  print(paste(vcount(g_used), ecount(g_used)))
  rank_networks_PPI[[i]] <- process_rank_network(network_set = "ppi", 
                                                 all_diseases = all_diseases,
                                                 el_all = list(ppi =el_used), 
                                                 weighted = F, disease_specific = F, to_exclude = NULL)
}


# Process fold results ###########################

rank_networks_processed <- lapply(rank_networks_PPI, function(x) 
  lapply(x, function(df) 
    lapply(1:10, function(iter) 
      df %>% filter(fold == iter) %>% select(trueset, rank))))


## compute AUC from the rank through CVAUC - combining cross validation results
pacman::p_load(cvAUC)
aucvals <- pbapply::pblapply(rank_networks_processed,  function(x) lapply(x, auc_output))


######
AUC_folds_each <- lapply(1:length(aucvals), function(x) auc_to_df(aucvals[[x]], label = names(aucvals)[x]))

AUC_folds <- bind_rows(AUC_folds_each)

saveRDS(AUC_folds, "../cache/fold_cv_processed_results_revision_LargeScalePPIRandomisation.RDS")

##################

# PPI large-scale

rank_networks_largescalePPI <- list()

el_PPI_largescale_full = process_edgelist(read_tsv("../data/network_edgelist_additional/PPI_LargeScale.tsv", col_names = c("A","B"), skip = 1, col_types = 'cc'))

rank_networks_largescalePPI[["el_PPI_largescale_full"]] <- process_rank_network(network_set = "ppi", 
                                                            all_diseases = all_diseases,
                                                            el_all = list(ppi =el_PPI_largescale_full), 
                                                            weighted = F, disease_specific = F, to_exclude = NULL)



rank_networks_largescalePPI_processed <- lapply(rank_networks_largescalePPI, function(x) 
  lapply(x, function(df) 
    lapply(1:10, function(iter) 
      df %>% filter(fold == iter) %>% select(trueset, rank))))

aucvals_largescalePPI <- pbapply::pblapply(rank_networks_largescalePPI_processed,  function(x) lapply(x, auc_output))


######
AUC_folds_largescalePPI <- auc_to_df(aucvals_largescalePPI[[1]], label = names(aucvals_largescalePPI)[1])

saveRDS(AUC_folds_largescalePPI, "../cache/fold_cv_processed_results_revision_largescalePPI_full.RDS")
