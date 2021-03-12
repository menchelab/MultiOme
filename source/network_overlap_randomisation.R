## randomised networks
## IzeBuphamalai, 11.03.2021

## Goals
### to get a reference distribution for networks, and their degree of overlap

# Formulation
#  For a list of grpah with V and E, from maximal poolcompute pairwise overlap


# 0 - set up -------------
source("../functions/readdata_functions.R")
library(igraph)


# 1 - load required files -----------
g = process_graph_data("../data/network_edgelists/")

all_nodes = sapply(g, function(x) V(x)$name) %>% unlist %>% unique 

A = matrix()

prob = 0.4
A = matrix(sapply(1:16, function(x) rbinom(1, size = 1, prob = 0.2)), nrow = 4, dimnames = list(letters[1:4], letters[1:4]))
B = matrix(sapply(1:16, function(x) rbinom(1, size = 1, prob = 0.2)), nrow = 4, dimnames = list(letters[1:4], letters[1:4]))

B_r = B
dimnames(B_r)[[1]] = sample(dimnames(B_r)[[1]], nrow(B_r))
dimnames(B_r)[[2]] = dimnames(B_r)[[1]]

gB_r = graph_from_adjacency_matrix(B_r)
g_A = graph_from_adjacency_matrix(A)
g_B = graph_from_adjacency_matrix(B)


g_A = g$`co-essential`
g_B = g$GOBP
g_C = g$GOMF



network_sim_df
###############
gs <- g
# 
# all_nodes <- lapply(gs, function(x) V(x)$name) %>% as_vector() %>% unique 
# ids = 1:length(all_nodes); names(ids) = all_nodes

# for(i in 1:length(gs)){
#   V(gs[[i]])$name <-  ids[V(gs[[i]])$name]
# }

ecounts <- sapply(gs, ecount)

node_index_randomisation <- function(graph){
  #V(graph)$name <- sample(names(V(graph)), vcount(graph), replace = F)
  name <- sample(names(V(graph)), vcount(graph), replace = F)
  #return(graph)
  return(name)
  
}

# maximum number of randomization
max_rand <- 30
#g_rand <- lapply(gs, function(graph) lapply(1:max_rand, function(x) node_index_randomisation(graph)))
name_rand <- lapply(gs, function(graph) lapply(1:max_rand, function(x) node_index_randomisation(graph)))


blank_result_mat <- matrix(NA, nrow = max_rand, ncol = max_rand)


overlap_randomisation <- function(g, name_gA, name_gB){
  g_A <- g[[name_gA]]
  g_B <- g[[name_gB]]
  
  union_mat <- blank_result_mat
  intersect_mat <- blank_result_mat
  
  for(i in 1:30){
    V(g_A)$name <- name_rand[[name_gA]][[i]]
    
    for(j in 1:30){
      V(g_B)$name <- name_rand[[name_gB]][[j]]
      
      intersect_mat[i,j] <- graph.intersection(g_A, g_B) %>% ecount
      union_mat[i,j] <- graph.union(g_A, g_B) %>% ecount
    }
  }
 
  return(list(intersect = intersect_mat %>% as.vector, union = union_mat %>% as.vector))
}


  
# graph_overlap <- function(g_A, g_B){
#   intersect <- apply(index_mat, 1, function(x) ) 
#   union  <- apply(index_mat, 1, function(x) graph.union(g_rand[[g_A]][[x[1]]], g_rand[[g_B]][[x[2]]]) %>% ecount ) 
#   return(list(intersect, union))
# }


pacman::p_load(pbapply)
network_sim_df$min_ecount <- pbapply(network_sim_df, 1, function(x) min(ecounts[x[1]], ecounts[x[2]]))


# this is the index matrix used for internally computing 
#ndex_mat <- tidyr::crossing(i = 1:max_rand, j = 1:max_rand)

randomisation_overlap_result <- pbapply(network_sim_df, 1, function(x) overlap_randomisation(g, x[1], x[2])) 

jaccard_index <- lapply(randomisation_overlap, function(x) x[[1]]/x[[2]])
overlap_index <- lapply(1:length(randomisation_overlap), function(x) randomisation_overlap[[x]][[1]]/graph_combn_df$min_ecount[x])
