### function to process passing layers of interest into propagation -------
source("../functions/readdata_functions.R")
source("../functions/LCC_functions.R")
source("../functions/process_edgelist.R")
source("../functions/process_LCC_result.R")
source("../functions/weighted_multiplex_propagation.R")
source("../functions/RWR_get_allnodes.R")
source("../functions/RWR_transitional_matrix.R")
source("../functions//RWR.R")

process_rank_network <- function(network_weight_df_custom = NULL, network_set=NULL, all_diseases, weighted = F, el_all, disease_specific = F, to_exclude = NULL){
  
  el = NULL 
  if(!disease_specific){
    # if a aprticular network, all a set of networks are used, unspecific to the disease group, el and supraadj are computed before passing into disease loop.
    
    if(is.null(network_weight_df_custom)){
      LCC_val_signif = all_LCC_val_signif %>%
        dplyr::filter(name %in% all_diseases, network %in% network_set, !network %in% to_exclude) %>%
        distinct(network)
      
      LCC_val_signif$LCC.zscore = 1
    } else{
      LCC_val_signif <- network_weight_df_custom
    }
    
    el = el_all[LCC_val_signif$network]
    
    gene_allnet = get_allnodes(el)
    supraadj = supraadjacency_compute(LCC_val_signif, el = el)
  }
  
  
  rank_df_all_folds = list()
  
  for(disease_current in all_diseases){
    
    print(disease_current)
    
    if(disease_specific){
      # if only consider significant network only, disease specific argument required
      
      LCC_val_signif = all_LCC_val_signif %>% 
        dplyr::filter(name == disease_current,  !network %in% to_exclude) %>% 
        select(network, LCC.zscore) 
      
      if(!weighted){LCC_val_signif$LCC.zscore = 1}
      
      print(sprintf("Number of network taken into account is %i", nrow(LCC_val_signif)))
      
      el = el_all[LCC_val_signif$network]
      gene_allnet = get_allnodes(el)
      supraadj = supraadjacency_compute(LCC_val_signif, el = el)
      
    }
    
    
    # set seed
    disease_genes = intersect(Orphanet_df[[disease_current]], gene_allnet)
    n = length(disease_genes)
    
    ## shuffle data
    seed = 144
    disease_genes = disease_genes[sample(n)]
    
    ## k-fold validation: k = 10
    k=10
    
    # subsets the disease genes for seed and retrieval sets
    set_labels = cut(1:n, breaks = k, labels = 1:k)
    
    
    rank_df = list()
    for(i in 1:k){
     # print(i)
      seed_sets = disease_genes[set_labels!=i]
      retrieval_sets = disease_genes[set_labels==i]
      ## perform retrieval
      
      ranks = weighted_multiplex_propagation(
        seedset = seed_sets, 
        trueset = disease_genes, 
        weighted_layer_df = LCC_val_signif,
        network_dirs = network_dirs, 
        el = el, 
        gene_allnet = gene_allnet, 
        Stest = supraadj)
      
      rank_df[[i]] = ranks[, colnames(ranks) %in% c("trueset", "rank", "RankArithmP")]
      colnames(rank_df[[i]]) <- c("trueset", "rank")
    }
    
    all_ranks<-  bind_rows(rank_df, .id = "fold") %>% mutate(max = length(gene_allnet))
    rank_df_all_folds[[disease_current]] = all_ranks
  }
  
  return(rank_df_all_folds)
}
