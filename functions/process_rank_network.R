### function to process passing layers of interest into propagation -------
source("../functions/readdata_functions.R")
source("../functions/LCC_functions.R")
source("../functions/process_edgelist.R")
source("../functions/process_LCC_result.R")
source("../functions/weighted_multiplex_propagation.R")
source("../functions/RWR_get_allnodes.R")
source("../functions/RWR_transitional_matrix.R")
source("../functions//RWR.R")

process_rank_network <- function(network_weight_df_custom = NULL, network_set=NULL, all_diseases, weighted = F, el_all, computeOnlyIfSignif = T, disease_specific = F, to_exclude = NULL){
  #' @computeOnlyIfSignif: the for a given disease, the network will be included if it was detected to be significant from the prior matric (LCC_vall_signif), if FALSE, this will be overridden and the retrieval performance will be computed for the disease regardless of whether the network was significant
  #if(computeOnlyIfSignif){
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
 # } #else{
    #el = network_specific_prop[["el"]]
    #gene_allnet = get_allnodes(el)
    #LCC_val_signif <- tibble(network = names(el), LCC.zscore = 1)
    #supraadj = supraadjacency_compute(LCC_val_signif, el = el)
  #}
  
  
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
    print(n)
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


############
# Process the AUC values
# function to compute AUC from rank
auc_output<-function(rank_df){
  #' @input rank_df: a two column data frame, with label and prediction respectively
  rank_label<-lapply(rank_df, function(x) x[,1])
  
  # check if 
  if( all(sapply(rank_label, function(x) any(x)))){
    rank_prediction<-lapply(rank_df, function(x) -x[,2])
    
    out <- cvAUC(rank_prediction, rank_label)
    return(out)
  } 
}


# process auc data into data frame for post analysis
# AUC plot for all diseases
auc_to_df <- function(auc_object, label = NULL){
  AUC_folds<-lapply(auc_object[!sapply(auc_object, is.null)], function(x) x$fold.AUC)
  
  AUC_folds<- reshape2::melt(AUC_folds)
  colnames(AUC_folds)<-c("AUC", "name")
  
  AUC_folds<-AUC_folds %>% mutate(name=fct_reorder(name, AUC))
  
  if(!is.null(label)){
    AUC_folds$label <- label
  }
  return(AUC_folds)
}

