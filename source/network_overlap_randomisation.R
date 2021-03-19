## network similarity randomisation
## IzeBuphamalai, 11.03.2021

## Goals
### to get a reference distribution of overlaps for all pairs of networks, and compute the significance

if(!file.exists("../cache/overlap_similarity_randomisation_results.RDS")){
  
  print("no precomputed results found. Starting the randomisation. This will take long...")
  
  # 0 - set up -------------
  library(pacman)
  
  p_load(reshape2, pbpply, tidyverse, igraph)
  
  source("../functions/readdata_functions.R")
  
  # load graphs
  g = process_graph_data("../data/network_edgelists/")
  
  
  # 1 - set up functions and preliminary variables -----------
  
  # a vector of all nodes on all layers
  all_nodes <- lapply(g, function(x) V(x)$name) %>% unlist %>% unique()
  
  # counts edges for each layers
  ecounts <- sapply(g, ecount)
  
  
  node_index_randomisation <- function(graph){
    # a function to randomise node ids
    # input: a graph
    # output: a randomised node ids
    
    name <- sample(names(V(graph)), vcount(graph), replace = F)
    
    return(name)
  }
  
  # set for number of network iterations (total pairs = max_rand^2)
  max_rand <- 10
  
  # a list storing randmomised node ids, each network with the legth of max_rand
  name_rand <- lapply(g, function(graph) lapply(1:max_rand, function(x) node_index_randomisation(graph)))
  
  
  # create a global variable for storing blank result for each iteration pair to be fed in the randomisation
  # matrix is blank and is of nxn, n = max_rand
  blank_result_mat <- matrix(NA, nrow = max_rand, ncol = max_rand)
  
  
  overlap_randomisation <- function(g, name_gA, name_gB, iter = max_rand){
    # a function to create a pair of randomised networks, and measure overlap and union between them: this function has internal for-loop, which was intentially made, to avoid multiple loading of redundant data
    # input: g: original graph list, name_gA/B: name index of the graph list, iter: number of iterations 
    # output: a list of intersect and union of edges from the randomised networks
    
    
    # create random graphs from indices A and B
    g_A <- g[[name_gA]]
    g_B <- g[[name_gB]]
    
    # create union and intersect matrices to store results from each computation
    union_mat <- blank_result_mat
    intersect_mat <- blank_result_mat
    
    # use for loops instead of sequential apply framework to save memories and data loading: e
    for(i in 1:iter){
      V(g_A)$name <- name_rand[[name_gA]][[i]]
      
      for(j in 1:iter){
        V(g_B)$name <- name_rand[[name_gB]][[j]]
        
        # compute intersect and union of graphs
        intersect_mat[i,j] <- graph.intersection(g_A, g_B) %>% ecount
        union_mat[i,j] <- graph.union(g_A, g_B) %>% ecount
      }
    }
    
    return(list(intersect = intersect_mat %>% as.vector, union = union_mat %>% as.vector))
  }
  
  
  
  # 2 - perform computation -----------
  network_sim_df <- readRDS("../cache/network_jaccard_overlap_similarity_df.RDS") %>% dplyr::filter(!grepl("core", V1), !grepl("core" , V2))
  
  
  randomisation_overlap_result <- pbapply(network_sim_df, 1, function(x) overlap_randomisation(g, x[1], x[2])) 
  
  # save the result in RDS format and to cache
  saveRDS(randomisation_overlap_result, "../cache/overlap_similarity_randomisation_results.RDS")
  
} else{
  
  print("load the precomputed value from cache")
  
  randomisation_overlap_result <- readRDS("../cache/overlap_similarity_randomisation_results.RDS")
  
}

