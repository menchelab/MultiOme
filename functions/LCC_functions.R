#' NetSig - measuring network signature of a node set of interest
#' Ize Buphamalai, CeMM, Vienna
#' Last updated: 18 Feb 2019

################################
#' LCC extract: create a function to extract nodes of interest out of the network and get the component size
#' 
#' @param g igraph object
#' @param node_set a vector of nodes of interest
#' @return ALL component sizes of node sets in the graph

lcc_subtract = function(g, gene_set){
  g_exclude = g - setdiff(V(g)$name, gene_set)
  lcc = components(g_exclude)$csize
  return(lcc)
}

# return igraph of the components
## save the igraph object of each subnetwork: layout=nicely
LCC_igraph_compute = function(graph, nodesets){
  g_CC = graph - setdiff(V(graph)$name, nodesets)
  return(g_CC)
} 

#  if(ecount(g_CC)>0){
#    l <- layout_nicely(g_CC)
#    sig <- sigmaFromIgraph(g_CC, layout = l)
#  }
#  else{
#    sig = NULL
#  }
#  return(sig)
  


################################
#' Measure LCC of a graph and perform randomisation control
#' 
#' @param graph an igraph object
#' @param nodesets node set to measure LCC 
#' @param trial the number of randomisation, default to 1000 
#' @param minnode the number of minimal number of node required for LCC to be performed. Otherwise is ignored
#' @return a vector containing:
#' (1) N: size of the node set in graph
#' (2) LCC: the size of the largest connected component
#' (3) mean, sd of the randomised gene set
#' (4) z - z-score computed to measure significance

LCC_randomisation_measure = function(graph, nodesets, minnode = 10, randomise = T, trial = 1000, node.randomisation = T){
  nodes = V(graph)$name
  nodesets_in_graph = intersect(nodes, nodesets)
  N = length(nodesets_in_graph)
  
  # there is no need to compute LCC if N < threshold, the result will not be taken into account in any case
  if(N>=minnode){
    LCC = max(lcc_subtract(graph, nodesets))
  }
  else{
    LCC = NA
  }
  
 if(randomise){
   # perform node or edge randomisation

   if(!is.na(LCC)){
     # perform the randomisation only if LCC was computed
     if(node.randomisation){ 
       # Performing node randomisation
       # rand_result = sapply(1:trial, function(k) sapply(N, function(x) max(lcc_subtract(graph, sample(nodes, x)))))
       rand_result = sapply(1:trial, function(k)  max(lcc_subtract(graph, sample(nodes, N))))
     } else{
       # use degree-preserving randomisation. Warning: slow and need a lot of iteration to properly work
       rand_result = sapply(1:trial, function(k)  max(lcc_subtract(rewire(graph, keeping_degseq()), nodesets_in_graph)))
     }
     
     mean = mean(rand_result)
     sd = sd(rand_result)
     z = (LCC-mean)/sd
   }
   
   else{
     # if LCC wasn't computed, the stat values were assigned as NA
     mean = NA
     sd = NA
     z = NA
   }
   
   
   
   return(c(N_in_graph = N, LCC.size = LCC, LCC.mean = mean, LCC.sd =sd, LCC.zscore =z))
 }
  else{
    return(c(N_in_graph = N, LCC.size = LCC))
  }
}

##################################
#' compute connectivity p-value as hypergeometric p-values of connectivity of all nodes with the seed in the graph
#' @param g an igraph object
#' @param node node to compute connectivity p-value
#' @param seed the set of seeds 
#' @return a vector containing:
# significantly connected to seed pval
connectvity_pval = function(g, node, seed){
  allnodes = V(g)$name
  if(!node %in% allnodes){
    vals = c(0,0,1)
    names(vals)= c("seed_neighbours", "all_neighbours", "pval")
  } else{
    node_neighbour = neighbors(g, node)$name
    N_seed = sum(seed %in% allnodes)
    N_neighbour = length(node_neighbour)
    N_seed_neighbour = sum(node_neighbour %in% seed)
    pval = phyper(N_seed_neighbour, N_seed, length(allnodes) - N_seed, N_neighbour, lower.tail = F)
    vals = c(N_seed_neighbour, N_neighbour, pval)
    names(vals)= c("seed_neighbours", "all_neighbours", "pval")
  }
  return(vals)
}
