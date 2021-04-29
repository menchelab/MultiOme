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

## from LCC_analysis_new: process the LCC
result_df = process_LCC_result(result_df)

### extract the most significant networks for all diseases
all_LCC_val_signif = result_df %>% dplyr::filter(correctedPval < 0.05) %>% dplyr::select(name, network, LCC.zscore)

###############################################
# read all edgelists for all significant networks
network_dirs = "../data/network_edgelists/"

el_all = list()
for(i in unique(result_df$network)){
  el_all[[i]] = read_tsv(paste0(network_dirs, i,".tsv"), col_names = c("A","B"), skip = 1)
}

el_all = lapply(el_all, process_edgelist)


############## 
# laod patients HPO terms
patient_hpo_df <- read_csv("../data/variants_list/patient_hpo_terms.csv", col_names = c('Patient', "hpo_id", "hpo_label"), skip = 1) %>%
  mutate(hpo_id = str_trim(hpo_id, side = 'both'))


if(!file.exists("../cache/patient_hpo_to_genes.RDS")){
  # create a function for retrieving genes associated to a hpo term
  hpo_genes_query = function(term){
    url_hpo <- httr::GET(paste0("https://hpo.jax.org/api/hpo/term/",term,"/genes?max=-1"))
    gene_list <- httr::content(url_hpo)
    
    # return gene label 
    gene_names <- sapply(gene_list$genes, function(x) x$entrezGeneSymbol)
    return(list(gene_names))
  }
  
  hpo_gene_association <- sapply(unique(patient_hpo_df$hpo_id), hpo_genes_query)
  
  saveRDS(hpo_gene_association, "../cache/patient_hpo_to_genes.RDS")
} else{
  hpo_gene_association <- readRDS("../cache/patient_hpo_to_genes.RDS")
}


# take only genes exist in network as seeds
seed_genes <- patient_hpo_df %>% 
  rowwise() %>%
  mutate(associated_genes = hpo_gene_association[hpo_id]) %>%
  summarise(Patient, genes = unlist(associated_genes))


#####
# read patient gene
gene_list <- read_csv("../data/patient_gene_list.csv") 



patient_annotated_df <- list()

# Perform prediction in three scenarios
# 1 all networks, 2 only significant networks to the disease group, 3 only ppi 


for(network_set in c("all", "signif", "ppi")){
  
  if(network_set == "signif"){
    
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
  for(patient_current in unique(patient_hpo_df$Patient)){
    print(patient_current)
    
    # select genes and weights for patient
    seeds <- seed_genes %>% 
      filter(Patient==patient_current, genes %in% gene_allnet) %>% 
      count(genes)
    
    
    ## perform retrieval
    ranks[[patient_current]] <- weighted_multiplex_propagation(
      seedset = seeds$genes, 
      seedweight = seeds$n,     
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
  
  patient_annotated_df[[network_set]] <- gene_list %>% 
    #dplyr::filter(Class=="H") %>%
    left_join(., ranks_all_patient_df[,c("Patient","GeneName","avg")], by = c("Patient", "GENE"="GeneName"))    %>%
    group_by(Patient) %>%
    distinct(GENE, .keep_all = T) %>%
    arrange(-avg) %>%
    mutate(network_rank = rank(-avg), total_variants = n())
}

###################
# merge results from all methods
patient_annotated_df_all_methods <- bind_rows(patient_annotated_df, .id = "network_set")

saveRDS(patient_annotated_df_all_methods, "../cache/patient_ranking_results.RDS")
