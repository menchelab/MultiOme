supraadjacency_compute = function(weighted_layer_df, el = NULL){
  if(is.null(el)){
    print("creating edge list from significant networks")
    el = list()
    for(i in weighted_layer_df$network){
      el[[i]] = read_tsv(network_dirs[grepl(i, network_dirs)], col_names = c("A","B"), skip = 1)
    }
    
    el = lapply(el, process_edgelist)
  }
  
  ############
  # Prepare the transitional matrix of all selected networks
  
  ## Procedures
  # 1. compute the supraadjajency - in sparse matrix form
  # 2. get the weights from network importances for the respective disease group (given)
  # 3. compare the performance with monolayer/ all inclusive layer without weights
  #############
  
  # transitional matrix for each layer
  library(Matrix)
 # print("creating supra transitional matrix")
  
  gene_allnet = get_allnodes(el)
  NLayer = length(el)
  namesLayer = names(el)
  
  # transitional matrix
  Mlist = lapply(el, function(x) transitional_matrix(x, gene_allnet))

  ## pmat: either taken from score obtained in network measurement or assign weights on them directly
  ## order similar to order in the el list in the beginning
  pmat = pmat_cal(weighted_layer_df$LCC.zscore)
  
  # supra-transitional matrix
  Stest = supratransitional(Mlist, pmat)
  
  return(Stest)
}

weighted_multiplex_propagation = function(seedset, seedweight = NULL, trueset, weighted_layer_df, network_dirs = NULL, el = NULL, gene_allnet = NULL, Stest = NULL, remove_seeds = TRUE){
  ###########
  #

    ##########
  require(igraph)
  require(tidyverse)

  
  if(is.null(network_dirs)){
    network_dirs = c(list.files("./networks/network_edgelists/",full.names = T)
    )
  }
  
  # report numbers from the input
  #print(sprintf("Number of seed used: %i, true set: %i", length(seedset), length(trueset)))

 # for(i in 1:nrow(weighted_layer_df)){
  #  print(paste0(" network layer: ", i, ": ", weighted_layer_df[i,1]))
  #}
  
  
  if(is.null(gene_allnet)){
    gene_allnet = get_allnodes(el)
  }
  
  
  NLayer = nrow(weighted_layer_df)
  namesLayer = weighted_layer_df$network
  
  if(is.null(Stest)){
    Stest = supraadjacency_compute(LCC_val_signif, el = el)
  }
  
  #############
  # add seed to the computation
  #############
  
  #print("Initialising the propagation")
  
  # create initial visiting probability: from seed
  p_0 = rep(0, length(gene_allnet))
  names(p_0) = gene_allnet
  
  # position of seed
  seed <- gene_allnet %in% seedset
  
  if(is.null(seedweight)){
    p_0[seedset] = 1
  } else{
    p_0[seedset] = seedweight
  }
  
  #detach names for memory efficiency
  names(p_0) = NULL
  
  p_0 = rep(p_0, NLayer)
  p_0 = p_0/sum(p_0)
  
    
  pos_nonseed = seq_along(gene_allnet)[!seed]
  
 # print("Performing the network-based propagation")
  test_result = RWR(M = Stest,p_0 = p_0, r = 0.7)
  result_df = matrix(test_result, ncol = NLayer, byrow = F,
                     dimnames = list(gene_allnet, namesLayer)) %>% as.data.frame()
  
  
  # annotate the data
  result_df$GeneName = rownames(result_df)
  result_df$seed = seed
  result_df$trueset = result_df$GeneName %in% as.character(trueset)
  
  
  #################
  # Combine rank values
  #################
  #summary stat if taken more than one layer into account
  if(NLayer > 1){
    ## 1. Arithmetic  mean on probability
    ######################################
    result_df$avg = rowMeans(result_df[,1:NLayer])
    
    ## 2. Geometric mean on probability
    ######################################
    result_df$GeomAvg = apply(result_df[,1:NLayer], 1, function(x) prod(x)^(1/NLayer)) 
    
    
    # remove seed from the result list
    if(remove_seeds){
      result_df_noseed = result_df %>% dplyr::filter(!seed)
    } else{
      result_df_noseed = result_df
    }
   
    
    # rank results
    
    ## 3. Rank all the nodes together (gene from different layers as different nodes and average their ranks)
    ######################################
    rank_alltogether_df = matrix(rank(-result_df_noseed[,1:NLayer]), ncol = NLayer, byrow = F)
    colnames(rank_alltogether_df) = namesLayer
    rank_alltogether_df = as.data.frame(rank_alltogether_df)
    rank_alltogether_df$GeneName = result_df_noseed$GeneName
    rank_alltogether_df$RankAllAvg = apply(rank_alltogether_df[,1:NLayer], 1, function(x) prod(x)^(1/NLayer)) 
    
    ## 4. Rank each layer separately and merge them together
    ######################################
    rank_eachlayer_df = apply(-result_df_noseed[,1:NLayer], 2, rank)
    rank_eachlayer_df = as.data.frame(rank_eachlayer_df) %>% mutate(GeneName = result_df_noseed$GeneName)
    rank_eachlayer_df$RankEachAvg = apply(rank_eachlayer_df[,1:NLayer], 1, function(x) prod(x)^(1/NLayer)) 
    
    #### combine rank from different sources
    ######################################
    
    rank_df = cbind(result_df_noseed[,-c(1:NLayer)], 
                    RankAllAvg= rank_alltogether_df$RankAllAvg ,
                    RankEachAvg =rank_eachlayer_df$RankEachAvg )
    
    rank_df = rank_df %>% mutate(RankArithmP = rank(-avg),
                                 RankGeomP = rank(-GeomAvg),
                                 RankAllAvg = rank(RankAllAvg),
                                 RankEachAvg = rank(RankEachAvg)
    )
    
  } else{
   # if there is only one layer, there is no summary statistics
    
    if(remove_seeds){
      result_df = result_df %>% dplyr::filter(!seed)  
      } 
    
    result_df$rank = rank(-result_df[,1])
    rank_df = result_df
  }

  # label if it gets to top 20 of any methods
  #rank_df$topgene = apply(rank_df[,6:9]<=50, 1, any)
  #rank_df$label = ifelse(rank_df$topgene, rank_df$GeneName, "")
  return(rank_df)
}
