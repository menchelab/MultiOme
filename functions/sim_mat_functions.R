# functions for network contruction

# remove n smallest values in the similarity matrix
# remove genes with minimum relevance - 
# - those with no connection at all and those that only connects via root 
# - i.e. first and second lowest value
# then clean up the matrix by remove genes with no connections

remove_minval = function(data, n){
  #function to remove n minimal value out of the original data framr
  for(i in 1:n){
    minval= min(data, na.rm = T)
    minval_roundUp = ceiling(minval*100)/100 #round up the second decimal place
    #cell_minval = data == minval
    cell_minval = data <= minval_roundUp
    message(sprintf("Iter %i, For weight %f, minimal relevance interaction removed: %i",
                    i, minval, sum(cell_minval, na.rm = T)))
    data[cell_minval] <- NA
  }
  
  # network has no self loop, so removing the diagonal items - replace by NA
  for(i in 1:nrow(data)){data[i,i] <- NA}
  ##
  message("remove those genes without any connections")
  k = apply(data, 2, function(x) sum(!is.na(x)))
  message(paste0(sum(k==0), " nodes removed")) # print out number of nodes removed
  # removing those nodes by rows and then by columns
  data = data[, k>0]
  data = data[k>0, ]
  return(data)
}

#######################################

# Backbone extraction
# input: value, p, k

#pval calculation for backbone - no need for actual integral

# backbone_cal = function(p, k){
#   if(is.na(p)){return(NA)} else{
#     int = integrate(function(x)(1-x)^(k-2), 0, p)
#     alpha = 1- (k-1)*int$val
#     return(alpha)
#   }
# }


# 1, getting matrix P, 
# where p_ij =w_ij/si = transitionl probability matrix = each cell divide by colsum
backbone_cal = function(weighted_adj_mat){
 pval_mat = weighted_adj_mat
 #====
 message("Calculate weight and degree of all nodes")
 W = colSums(weighted_adj_mat, na.rm = T)
 k = apply(weighted_adj_mat, 2, function(x) sum(!is.na(x)))
 #====
 message("Calculate p-value from each connection")
 for(i in 1:ncol(pval_mat)){
    pval_mat[,i] = (1-(weighted_adj_mat[,i]/W[i]))^(k[i]-1)
 }
 return(pval_mat)
}

# Create edge list from that p value matrix
adj_mat_from_pvalmat = function(pvalmat, threshold){
  result_mat = pvalmat < threshold
  # adjacency matrix, obtained when both p[i.j] and p[j,i] match the criteria
  new_adj = t(result_mat)*result_mat
  new_adj[is.na(new_adj)] = 0
  # remove genes with no significance
  k = apply(new_adj, 2, sum)
  new_adj = new_adj[,k>0]
  new_adj = new_adj[k>0,]
  
  return(new_adj)
}

# write a function to plot the degree distribution
# https://chengjunwang.com/web_data_analysis/demo2_simulate_networks/
# plot and fit the power law distribution
fit_power_law = function(graph, mode = "pdf", title = "Degree Distribution") {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  #CDF format
  ccd = degree.distribution(graph, mode = "all", cumulative = T)
  cd = 1-ccd #complementary CDF
  degree = 1:max(d)
  if(mode == "pdf"){
    probability = dd[-1]; ylab = "pdf(k)"
  }
  else if(mode == "cdf"){
    probability = cd[-1]; ylab = "cdf(k)"
  }
  else if(mode == "ccdf"){ #complementary CDF
    probability = ccd[-1]; ylab = "ccdf(k)"
  }
  # delete blank values
  nonzero.position = which(dd[-1] != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  # plot
  plot(probability ~ degree, log = "xy", xlab = "k", ylab = ylab, 
       main = title, pch = 16, col = "grey50",
       sub = paste("Alpha =", round(alpha, 3), "; R sq =", round(R.square, 3)))
  curve(power.law.fit, col = "red", add = T, n = length(d))
}


## vectorise_edge_list
vectorised_edges = function(g){
  # function produce edge list in string
  elist_string = function(x1,x2){ result = paste(sort(c(x1, x2)), collapse ="-") 
  return(result)}
  # create edge data frame from graph
  elist = as.data.frame(as_edgelist(g))
  # convert edge lists to string
  name = apply(elist, 1, function(x) elist_string(x[[1]], x[[2]]))
  return(name)
}