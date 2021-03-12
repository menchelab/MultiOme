## Network similarity analysis
## basic plots of nodes and edges
## Ize Buphamalai
## update Aug 25, 2020


# 0 - set up -------------
source("../functions/readdata_functions.R")
library(igraph)
library(tidyverse)

# 1 - load required files -----------
g = process_graph_data("../data/network_edgelists/")

# add the pan-tissue co-expression network - co-expression networks found in over 5 tissues
coex_raw_edge_count = readRDS("../cache/coexpression_raw_edge_counts.RDS")
g$coex_core = coex_raw_edge_count %>%
  dplyr::filter(n>=5) %>% 
  select(A,B) %>% 
  graph_from_data_frame(., directed = F)

# remove large object from memory
rm(coex_raw_edge_count)

# 2 - defines the graph similarity function  -----------
graph_similarity = function(g, name1, name2){
  # compute the Jaccard and overlap Index
  
  
  g_intersect = graph.intersection(g[[name1]], g[[name2]], keep.all.vertices = FALSE)
  g_union = graph.union(g[[name1]], g[[name2]], byname = T)
  
  #Jaccard index
  jAB = ecount(g_intersect)/ecount(g_union)
  
  # Overlap index
  oAB = ecount(g_intersect)/min(ecount(g[[name1]]), ecount(g[[name2]]))
  return(c(jaccardIndex = jAB, overlapindex = oAB))
}

# 3- compute jaccard and overlap similarity ---------

# pairwise combination of all networks names
network_sim_df = combn(names(g),m = 2, simplify = T) %>% t() %>% as_tibble

# compute pairwise similarity - this process might take really long
simIndex = apply(network_sim_df, 1, function(x) graph_similarity(g, x[1], x[2])) 
simIndex_df <- as_tibble(t(simIndex))

# merge the results
network_sim_df <- bind_cols(network_sim_df, simIndex_df)

saveRDS(network_sim_df, file = "../cache/network_jaccard_overlap_similarity_df.RDS")

