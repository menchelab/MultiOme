#### Multiplex RWR for all the Orphanet: using k-fold cross validation, actually recomputing all the networks
## Ize Buphamalai
## 16 Apr 2020
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


## read all edgelists for all significant networks
network_dirs = "../data/network_edgelists/"
g = process_graph_data("../data/network_edgelists/")

el_all = lapply(g, function(x) process_edgelist(as_data_frame(x)))


# Read the disease gene association
################
Orphanet_df = process_disease_genes_data("../data/table_disease_gene_assoc_orphanet_genetic.tsv", min_gene = 20, max_gene = 2000)
Orphanet_df = Orphanet_df$disgene_list

all_diseases = names(Orphanet_df)



# load the precomputed networks
lcc_stat_vals  = readRDS("../cache/LCC_stat_compute_vals_all_networks.RDS")
lcc_stat_vals_full_df = bind_rows(lcc_stat_vals, .id = "network") %>% dplyr::filter(type == "predicted") %>% select(-type)
#lcc_stat_vals_full_df = left_join(lcc_stat_vals_full_df, network_info)



## to store ranking results
rank_df_all_folds = list()

network_signif_all_folds = list()

#profvis({
  for(disease_current in all_diseases){
    
    print(disease_current)
    
    # all disease genes in the set
    disease_genes = Orphanet_df[[disease_current]]
    n = length(disease_genes)
    
  
    ## shuffle data
    seed = 144
    disease_genes = disease_genes[sample(n)]
    
    ## k-fold validation: k = 10
    k=10
    
    # subsets the disease genes for seed and retrieval sets
    set_labels = cut(1:n, breaks = k, labels = 1:k)
    
    
    # each retrievat sets controlled by parameter i, i = 1:k
    rank_df = list()
    signif_nets = list()
    
    for(i in 1:k){
      seed_sets = disease_genes[set_labels!=i]
      retrieval_sets = disease_genes[set_labels==i]
      
      # if the precomputed LCC exists, don't recompute
      if(!is.null(network_signif_all_folds[[disease_current]])){
        LCC_result = sapply(g, function(x) LCC_randomisation_measure(graph = x, nodesets = seed_sets, randomise = F))
        result_df = t(LCC_result) %>% as.data.frame() %>% rownames_to_column(., var = "network")
        
        # add the mean and sd
        result_df_pass = result_df %>% dplyr::filter(!is.na(LCC.size)) %>% 
          left_join(., lcc_stat_vals_full_df, by = c("network", "N_in_graph" = "N_node")) %>%
          mutate(LCC.zscore = (LCC.size - mean)/sd,
                 name = disease_current)
        
        
        result_df_pass = process_LCC_result(result_df_pass, network_annotate = F)
        
        # store the significant levels for further analytics
        signif_nets[[i]] = result_df_pass
      }
      else{
        result_df_pass = network_signif_all_folds[[disease_current]] %>% dplyr::filter(iter == i)
      }
    
      ### extract the most significant networks for all diseases
      LCC_val_signif = result_df_pass %>% 
        dplyr::filter(correctedPval < 0.05) %>%
        dplyr::select(network, LCC.zscore)
    
    
      ### build edge list
      
      el = el_all[LCC_val_signif$network]
      gene_allnet = get_allnodes(el)
      supraadj = supraadjacency_compute(LCC_val_signif, el = el)
      
      
      ## perform retrieval
      
      rank_df[[i]] = weighted_multiplex_propagation(
        seedset = seed_sets, 
        trueset = disease_genes, 
        weighted_layer_df = LCC_val_signif,
        network_dirs = network_dirs, 
        el = el, 
        gene_allnet = gene_allnet, 
        Stest = supraadj)
    }
    
    
    
    rank_df_all_folds[[disease_current]] = rank_df
    #network_signif_all_folds[[disease_current]] = bind_rows(signif_nets, .id = "iter")
  }
  
#})

saveRDS(rank_df_all_folds, "./cache/LCC_CV_ranking_rare_genetic_diseases.RDS")
saveRDS(network_signif_all_folds, "./cache/LCC_network_significance_CV.RDS")

########
#plot using cvAUC

pacman::p_load("cvAUC")

# store AUC computed values
aucvals = list()

pdf(paste0("./Figs/cvROC/all_rare_genetic_diseases.pdf"), width = 12, height = 10)
par(mfrow=c(5,6))
{
  for(disease_current in all_diseases){
    rank_label = lapply(rank_df_all_folds[[disease_current]], function(x) x$trueset)
    rank_prediction = lapply(rank_df_all_folds[[disease_current]], function(x) -x$RankGeomP)
    
    out <- cvAUC(rank_prediction, rank_label)
    
    aucvals[[disease_current]] = out
    
    #Plot fold AUCs
     { plot(out$perf, col="grey82", lty=3, main=paste(disease_current))
       
       #Plot CV AUC
       plot(out$perf, col="red", avg="vertical", add=TRUE)}
  }
}

dev.off()


