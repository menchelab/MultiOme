## Patient specific RWR
# Ize Buphamalai, CeMM
# edited: 3 April 2021

source("../functions/readdata_functions.R")
source("../functions/LCC_functions.R")
source("../functions/process_edgelist.R")
source("../functions/process_LCC_result.R")
source("../functions/weighted_multiplex_propagation.R")
source("../functions/RWR_get_allnodes.R")
source("../functions/RWR_transitional_matrix.R")
source("../functions//RWR.R")
library(igraph)


## load the pre-computed LCC results 
result_folder = "../cache/output/Orphageneset_rare/"
result_df = readRDS(paste0(result_folder, "LCC_and_distance_calculation_results.RDS"))

result_df_additional_net = readRDS("../cache/output/Orphageneset_rare_additionalNetworks/LCC_and_distance_calculation_results.RDS")

result_df <- rbind(result_df,result_df_additional_net %>% filter(network == "coex_core"))

## from LCC_analysis_new: process the LCC
result_df = process_LCC_result(result_df)

### extract the most significant networks for all diseases
all_LCC_val_signif = result_df %>% dplyr::filter(correctedPval < 0.05) %>% dplyr::select(name, network, LCC.zscore)

###############################################
# read all edgelists for all significant networks
network_dirs = "../data/network_edgelist_combn//"

el_all = list()
for(i in unique(result_df$network)){
  el_all[[i]] = read_tsv(paste0(network_dirs, i,".tsv"), col_names = c("A","B"), skip = 1, col_types = 'ss')
}

el_all = lapply(el_all, process_edgelist)

# 

patient_annotated_df <- list()

# Perform prediction in three scenarios
# 1 all networks, 2 only significant networks to the disease group, 3 only ppi 


for(network_set in c("all", "informedMultiPlex", "ppi")){
  print(network_set)
  
  if(network_set == "informedMultiPlex"){
    
    LCC_val_signif = all_LCC_val_signif %>% 
      dplyr::filter(name == "Rare genetic neurological disorder") %>% 
      select(network, LCC.zscore)
    
    el = el_all[LCC_val_signif$network]
    gene_allnet = get_allnodes(el)
    supraadj = supraadjacency_compute(LCC_val_signif, el = el)
    
  } else if(network_set == "all"){
    el = el_all
    gene_allnet = get_allnodes(el_all)
    LCC_val_signif = tibble(network = unique(result_df$network), LCC.zscore = 1)
    supraadj = supraadjacency_compute(LCC_val_signif , el = el)
    
  } else if(network_set == "ppi"){
    el = el_all["ppi"]
    gene_allnet = get_allnodes(el)
    LCC_val_signif = tibble(network = "ppi", LCC.zscore = 1)
    supraadj = supraadjacency_compute(  LCC_val_signif, el = el)
  }
  
  
  ranks = list()
  for(patient_current in unique(patient_hpo_df_all_cohort$`Experiment ID`)){
    print(patient_current)
    
    # select genes and weights for patient
    seeds <- unlist(seed_genes$associated_genes[seed_genes$`Experiment ID`==patient_current])
    seeds <- tibble(genes = seeds[seeds %in% gene_allnet]) %>% count(genes)
    
    if(network_set == "informedMultiPlex"){
      seed_weight <- seeds$n
    } else{
      seed_weight <- NULL
    }
    
    ## perform retrieval
    ranks[[patient_current]] <- weighted_multiplex_propagation(
      seedset = seeds$genes, 
      seedweight = seed_weight,
      trueset = c(), 
      weighted_layer_df = LCC_val_signif,
      network_dirs = network_dirs, 
      el = el, 
      gene_allnet = gene_allnet, 
      Stest = supraadj,
      remove_seeds = F)
  }
  
  ranks_all_patient_df <- bind_rows(ranks, .id = "Patient")
  
  if(network_set=="ppi"){
    ranks_all_patient_df$avg = ranks_all_patient_df$ppi
  }
 
  rownames(ranks_all_patient_df) <- NULL
  ranks_all_patient_df$Patient <- as_factor(ranks_all_patient_df$Patient)
  ranks_all_patient_df$GeneName <- as_factor(ranks_all_patient_df$GeneName)
 
  
  patient_annotated_df[[network_set]] <- ranks_all_patient_df
  
# this is to be obtaibed when having full variant list   
#  patient_annotated_df[[network_set]] <- gene_list %>% 
    #dplyr::filter(Class=="H") %>%
#    left_join(., ranks_all_patient_df[,c("Patient","GeneName","avg")], by = c("Patient", "GENE"="GeneName"))    %>%
#    group_by(Patient) %>%
#    distinct(GENE, .keep_all = T) %>%
#    arrange(-avg) %>%
#    mutate(network_rank = rank(-avg), total_variants = n())
}

saveRDS(patient_annotated_df, "../cache/test_rdconnect_ranking_results.RDS")



