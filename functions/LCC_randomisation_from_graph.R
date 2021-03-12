## function to compute mean and SD of a network node randomization
## Ize Buphamalai
## 25.11.19

LCC_stat_compute = function(g, N_iter = 1000){
  #' @input g: an igraph object
  #' @output LCC_vals: a data frame with three columns: N_nodes, mean and sd
  require(pacman)
  pacman::p_load(pbapply, igraph, tidyverse)
  
  nodes = V(g)$name
  
  # LCC rand function to perform 1000 iterations of LCC measurements of a node set of size n_node
  LCC_rand = function(n_node){
    rand_result = sapply(1:N_iter, function(k)  max(lcc_subtract(g, sample(nodes, n_node))))
    mean = mean(rand_result)
    sd = sd(rand_result)
    return(c(mean, sd))
  }
  
  # node size to actually compute the values from: due to fluctuation at early node size: 
  # up to N=50 sampled with increment of 1, incremaent of 5 for 50 < N < 200, and increment of 10 for 200 < N < 2000
  #sequences = c(5:49, seq(50, 195, 5), seq(200, 2000, 10))
  breaks = c(10,50,100,500,1000,2050)
  increments = c(2, 5, 10, 20, 50)
  sequences = sapply(1:length(increments), function(i) seq(breaks[i], breaks[i+1]-increments[i], increments[i]))
  sequences = unlist(sequences)
  rand_results = pbsapply(sequences, LCC_rand)
  
  rand_result_df = as.data.frame(t(rand_results)); colnames(rand_result_df) = c("mean", "sd")
  rand_result_df$N_node = sequences
  rand_result_df$type = "computed"

  # interpolate N_nodes that were not sampled for calculation
  # mean and sd were computed using loess (local regression) model in R, with span of 10% neighbouring points and fit polynomial degree 2
  points_to_predict = setdiff(10:2000, rand_result_df$N_node)
  points_to_predict = 10:2000
  
  mean.lo = loess(mean ~ N_node, rand_result_df, span = 0.1, degree = 2)
  pred.mean = predict(mean.lo, data.frame(N_node = points_to_predict))
  
  sd.lo = loess(sd ~ N_node, rand_result_df, span = 0.1,degree =2)
  pred.sd = predict(sd.lo, data.frame(N_node = points_to_predict))
  
  predicted_results = data.frame(mean = pred.mean, sd = pred.sd, N_node = points_to_predict, type = "predicted")
  
  # combine predicted and computed data
  comb = rbind(rand_result_df, predicted_results)    
  comb = comb[order(comb$N_node),]
  
  # round mean and sd to max 3 digits
  comb$mean = round(comb$mean, 3)
  comb$sd = round(comb$sd, 3)
  
  # rearrange columns: N-node, mean and SD, and type
  LCC_vals = comb[, c(3, 1,2, 4)]
  
  return(LCC_vals)
}
