################################
#' compute the closest distance for a set of nodes. Note that it excludes distance to self (e.g. d=0)
#' this function is required by the distance calculation main function
#' 
#' @param d distance matrix of nodes to want to measure distance from
#' @return a vector of closest distance
d_closest = function(d){
  # check if it is matrix 
  if(is.null(dim(d))){
    d = matrix(d)
  }
  # compute closest distance
  d_c= apply(d, 2, function(x) min(x, na.rm = T))
  d_c = d_c[is.finite(d_c)]
  return(d_c)
}


##############################
#' A function to calculation distance matrix from an igraph object, with removing isolated nodes
#' 
#' @param graph an igraph object
#' @return distance matrix
dismat_cal = function(graph, ignore.isolate = T, ignore.self = F){
  dismat = distances(graph)
  
  # remove self-distance (diagonal) values (in case of intraset measures)
  if(ignore.self){
    dismat[!is.finite(dismat)] = NA
  }
  
  # detect isolated nodes
  if(ignore.isolate){
    rm_nodes = which(apply(dismat, 1, function(x) all(is.na(x))))
    }
  
  # remove isolated nodes
  if(length(rm_nodes)!=0){
    dismat = dismat[-rm_nodes, -rm_nodes]
  }
  return(dismat)
}


############################
#' Calculating the significance of the average distance measurement through z-score and fold change
#' 
#' @param val the measured average distance value
#' @param random_vals a numerical array of average distance of random node sts
#' @return signif_measure: a list consisting of zscore, corresponding pval, ans the odd ratio
signif_distance_cal = function(val, random_vals){
  mean_rand = mean(random_vals)
  sd_rand   = sd(random_vals)
  
  # statistical measures
  signif_measures = list()
  signif_measures[["zscore"]] = (val - mean_rand)/sd_rand
  signif_measures[["zpval"]] = pnorm(signif_measures[["zscore"]], lower.tail = T)
  signif_measures[["foldchange"]] = mean_rand/val
  return(signif_measures)
}


################################
#' Average distance calculation among nodes in the gene set (given one node set provided) or between two sets (given two node sets provided)
#' The distance can be computed in two modes e.g. 'closest' or 'shortest'
#' 
#' @param distmat distance matrix of a graph
#' @param g an igraph object, used as alternative input to distance matrix in case it is not provided. In this case, distance matrix will be included
#' @param vset1 vertex set 1
#' @param vset2 (optional) vertex set 2. It is NULL by default, which means the distance is only measured among nodes in set 1.
#' @param mode 'closest' or 'shortest'or both. 
#' @return a list consisting of
#' (1) number of nodes that the distances were measured from
#' (2) raw values of all distance measured between the two node set
#' (3) mean value of the measured distance 
#' Otherwise the function will measure distances between two node sets
distance_cal = function(g, dismat=NULL, vset1, vset2 =NULL , mode = c("closest", "shortest"), silent = F){
  
  # check if the distance matrix is already provided
  if(is.null(dismat)){
    if(!silent) print("distance matrix not provided, calculating...")
    dismat = dismat_cal(g, ignore.isolate = T, ignore.self = F)
  }
  
  # extract node lists in the graph
  nodes_in_graph = dimnames(dismat)[[1]]
  
  results = list()
  
  
  vset1 = vset1[vset1 %in% nodes_in_graph]
  
  # if there is no node in graph, quit the calculation
  if(length(vset1)>0){
    if(is.null(vset2)){
      if(!silent) message("perform the distance calculation within the node set")
      # in this case, it is internal comparison of a gene set
      
      # store the dimension of the node set
      results[["dim"]] = c(length(vset1), length(vset1))
      
      diag(dismat) = NA # prevent the calculatiion of self-distance
      
      # extract the distance matrix of subset of node lists
      d = dismat[vset1, vset1]
      
      # closest distance calculation
      ##################
      if("closest" %in% mode){
        if(!silent) message("performing closest distance calculation")
        
        results[["d_c"]] = d_closest(d)
        results[["d_c_mean"]] = mean(results[["d_c"]])
        
        # random closest distance is obtained by sampling random nodes
        # default sampling size of 1000
        
        d_c_rand_mean = c()
        for(i in 1:1000){
          rand_vset1 = sample(nodes_in_graph, length(vset1)) 
          d_rand = dismat[rand_vset1, rand_vset1]
          d_c_rand = d_closest(d_rand)
          d_c_rand_mean = append(d_c_rand_mean, mean(d_c_rand))
        }
        
        results[["d_c_rand"]] = d_c_rand_mean
        
        # significance calculation
        signif = signif_distance_cal(results[["d_c_mean"]], results[["d_c_rand"]])
        
        # append the values to the results
        results = append(results, unlist(list(d_c_signif = signif)))
      }
      
      # shortest distance calculation
      ##################
      if("shortest" %in% mode){
        if(!silent) message("performing shortest distance calculation")
        # in calculating the shortest pairwise distance, the self-distance (diagonal) is automatically excluded
        results[["d_s"]] = d[upper.tri(d, diag = F) & is.finite(d)]
        results[["d_s_mean"]] = mean(results[["d_s"]])
        
        d_s_random = dismat[upper.tri(dismat, diag = F) & is.finite(dismat)]
        
        # significance calculation
        signif = signif_distance_cal(results[["d_s_mean"]], d_s_random)
        
        # append the values to the results
        results = append(results, unlist(list(d_s_signif = signif)))
        
      }
    }
    
    if(!is.null(vset2)){
      # in this case, distances are calculated between nodes in the two sets
      if(!silent) message("perform the distance calculation between two node sets")
      
      vset1 = vset1[vset1 %in% nodes_in_graph]
      vset2 = vset2[vset2 %in% nodes_in_graph]
      results[["dim"]] = c(length(vset1), length(vset2))
      
      if(length(vset2)>0){
        d = dismat[vset1, vset2]
        
        # if vset1 or two has only one gene, d will be a vector and need to be turned into matrix before the next step
        if(!is.matrix(d)){
          d = as.matrix(d)
        }
        
        # remove isolate nodes (no distance to other nodes)
        isolate_row = apply(d, 1, function(x) all(is.infinite(x)))
        isolate_col = apply(d, 2, function(x) all(is.infinite(x)))
        d = d[!isolate_row, !isolate_col]
        
        if(!is.null(dim(d)) & "closest" %in% mode){
          if(!silent) message("performing closest distance calculation")
          d_c_row = apply(d, 1, function(x) min(x[is.finite(x)]))
          d_c_col = apply(d, 2, function(x) min(x[is.finite(x)]))
          results[["d_c"]] = c(d_c_row, d_c_col)
          results[["d_c_mean"]] = mean(c(d_c_row, d_c_col))
          
          # random closest distance is obtained by sampling random nodes
          d_c_rand_mean = c()
          for(i in 1:1000){
            rand_vset1 = sample(nodes_in_graph, length(vset1)) 
            rand_vset2 = sample(nodes_in_graph, length(vset2)) 
            d_rand = dismat[rand_vset1, rand_vset2]
            
            # remove isolate nodes (no distance to other nodes)
            isolate_row = apply(d_rand, 1, function(x) all(is.infinite(x)))
            isolate_col = apply(d_rand, 2, function(x) all(is.infinite(x)))
            d_rand = d_rand[!isolate_row, !isolate_col]
            
            d_c_rand = d_closest(d_rand)
            d_c_rand_mean = append(d_c_rand_mean, mean(d_c_rand))
          }
          results[["d_c_rand"]] = d_c_rand_mean
          
          # significance calculation
          signif = signif_distance_cal(results[["d_c_mean"]], results[["d_c_rand"]])
          
          # append the values to the results
          results = append(results, unlist(list(d_c_signif = signif)))
        }
        
        if(!is.null(dim(d)) & "shortest" %in% mode){
          if(!silent) message("performing shortest distance calculation")
          results[["d_s"]] = d[is.finite(d)] #same nodes between two sets count
          results[["d_s_mean"]] = mean(results[["d_s"]])
          d_s_random = dismat[upper.tri(dismat, diag = F) & is.finite(dismat)]
          
          # significance calculation
          signif = signif_distance_cal(results[["d_s_mean"]], d_s_random)
          
          # append the values to the results
          results = append(results, unlist(list(d_s_signif = signif)))
        }
      }
    }
  }
  

  return(results)
}


################################
#' Load distance matrix (if precalculated and stored in cache), or compute the matrix and store in cache
#' @param graph_name name of graphs to load distance matrix from
#' @return a distance matrix (pairwise)
load_dismat = function(distance_matrix_path = "./cache/distance_matrices/dismat_", prefix = NULL, graph_name){
  # Pre compute distance matrix
  if(!is.null(prefix)){
    dismat_path = paste0(distance_matrix_path, prefix, "_", graph_name,".RDS")
  } else{
    dismat_path = paste0(distance_matrix_path, graph_name,".RDS")
  }
  
  # if no distance matrix file yet, compute one  
  if(!file.exists(dismat_path)){
    # pre-calculate distance matrix
    message(paste("Computing distance matrix for",graph_name))
    dismat = dismat_cal(g[[graph_name]])
    saveRDS(dismat, dismat_path)
  } else{
    # otherwise load the precomputed matrix
    dismat = readRDS(dismat_path)
    message(paste("Loading the pre-computed distance matrix for ", graph_name))
  }
  return(dismat)
}
