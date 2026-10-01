# Verbatim copies of the original MultiOme R functions (git commit cce867d), used only
# to generate golden values for the Python parity tests. Do not edit.
# Sources: functions/RWR_transitional_matrix.R, functions/RWR.R, functions/RWR_get_allnodes.R

#' get transitional matrix across all networks
#' @param el A list of edge lists for different networks.
#' @param allnodes a vector containing unique names of all networks
#' @return  a sparse transitional matrix
#' @examples
#' el_list = list(edgelist1, edgelist2, edgelist3)
#' allnodes= get_allnodes(el_list)
#' Mlist = lapply(el_list, function(x) transitional_matrix(x, gene_allnet))
transitional_matrix = function(el, allnodes){
  #' calculate the adjacency matrix based on edge list
  #' @input: a list of edgelist, allnodes
  
  # position of the elements to be one
  posA = as.numeric(factor(el$A, levels = allnodes))
  posB = as.numeric(factor(el$B, levels = allnodes))
  n = length(allnodes)
  # extends the position to the last element of the matrix- which will be set to zero
  #vals = c(rep(1, length(posA)), 0)
  vals = rep(1, length(posA))
  #posA = c(posA, n); posB = c(posB, n)
  
  A = sparseMatrix(i=posA, j=posB, x = vals, 
                   dimnames = list(allnodes, allnodes), symmetric = T, dims = c(n,n))
  M = A %*% Matrix::Diagonal(x = 1 / Matrix::colSums(A))
  return(M)
}


#' get inter-later transitional matrix based on weight vector assigned for each layer, computed based on the principle of detailed balance 
#' @param w A weight vector for each layer
#' @return  inter-layter transitional matrix
#' @examples
#' S = transitional_matrix(el_list, allnodes)
pmat_cal = function(w){
  L = length(w)
  p = matrix(0 , nrow = L, ncol = L)
  for(i in 1:L){
    for(j in setdiff(1:L,i)){
      p[i,j] = (min(1, w[i]/w[j]))/(L)
    }
  }
  diag(p) = 1 - colSums(p)
  return(p)
}


#' get supra transitional matrix based on single-layered transitional matrix
#' @param Mlist A list of transitional matrices for each layer
#' @param pmat a vector containing weights for each layer 
#' @return  a sparse supra-transitional matrix
#' @examples
#' S = transitional_matrix(el_list, allnodes)
supratransitional = function(Mlist, pmat){
  n_gene = nrow(Mlist[[1]])
  nL = length(Mlist)
  for(i in 1:nL){
    
    for(j in 1:nL){
      
      # determine the block i,j value if it is the intra- or interlayer transitional
      # interlayer transitional is simply the diagonal
      if(i == j){blockmat = Mlist[[i]]} else blockmat = Diagonal(n_gene) 
      # adjust the blockmat with the weight
      weighted_blockmat = pmat[i,j]*blockmat
      
      # the supra-transitional is adding up each row separately:
      # if j=1 (first column, the beginning of the block), assign it to the block matrix, otherwise append it to the right
      if(j==1){rowblock = weighted_blockmat}
      else rowblock = cbind(rowblock, weighted_blockmat)
    }
    
    # add the rowblock to the existing ones
    if(i==1){S = rowblock}
    else S = rbind(S, rowblock)
   # print(paste0(i, " out of ", nL, " layers completed"))
  }
  
  # correct for the case for nodes don't exist in all layers, leading to ColSums for some elements less than 1. 
  #Normalising it all would solve the problem
  S = S %*% Matrix::Diagonal(x = 1 / Matrix::colSums(S))
  return(S)
}



RWR <- function(M, p_0, r=0.8 , prop=FALSE, scale = F) {
  
  require("Matrix")

  # use Network propagation when prop=TRUE
  if(prop) {
    w<-Matrix::colSums(M)
    # Degree weighted adjacency matrix
    # define weight matrix, w(i,j) = w(j,i) = sqrt(ki*kj)
    
    W<-matrix(rep(1,nrow(M)*ncol(M)), nrow(M), ncol(M))
    for (i in 1:nrow(M)) {
      for (j in 1:ncol(M)) {
        v<-sqrt(w[i]*w[j])
        W[i,j]<-v
        W[j,i]<-v
      }
    }
    
    # normalise the adj matrix by weight
    W<-M/W
    
    # in the case of non propagation
  } else if(scale) {
    # Convert M to column-normalized adjacency matrix (transition matrix, depend only on the outgoing nodes)
    W=scale(M,center=F,scale=Matrix::colSums(M))
  }
  else{W = M}
  
  # Assign equal probabilities to seed nodes
  
  p_0<-p_0/sum(p_0)
  p_t<-p_0
  # Iterate till convergance is met
  converge = F
  D_t = c(1000, 100) # just some initial large values
  i=1
  while (!converge ) {
    # Calculate new proabalities
    p_tx <- (1-r) * W %*% p_t + r * p_0
    # Check convergance
    D_tx = norm(p_tx-p_t)
    D_t = append(D_t, D_tx)
    #print (D_tx)
   
     # convergence is taken if the values don't change, or change very minimally (ratio test)
    if(abs(D_t[i+2]-D_t[i+1])==0){
      converge = T
    } else {
       conv_ratio = log10(abs(D_t[i+2]-D_t[i+1]))/log10(abs(D_t[i+1]-D_t[i]))
    if ( round(conv_ratio, 2) %in% c(0.99,1,1.01) ) {
      converge = T
    } else {
      # converge if lopped for over a hundred times
      if(i >100)
        converge = T
    } #else 
    #f(abs(D_t[i+2]-D_t[i+1]) < 1e-23){
    #  converge = T
    #}
      
       }
    #else{
      #D_t = D_tx
    #}
    p_t<-p_tx
    i= i+1
  }
  
  return(p_tx)
}

#' get a vector of names of all nodes from all networks
#' 
#' @param el_list A list of edge lists for different networks.
#' @return  a vector containing unique names of all networks, required for building adjacency matrix
#' @examples
#' el_list = list(edgelist1, edgelist2, edgelist3)
#' allnodes= get_allnodes(el_list)
#' get_allnodes = function(el_list){
#'   #' get all nodes from several edge lists
#'   # merge edgelist 
#'   el_merged = do.call(rbind.data.frame, el_list)
#'   el_merged = el_merged %>% count(A,B)
#'   
#'   # gene_allnet = unique(c(unique(el_merged$A), unique(el_merged$B)))
#'   gene_allnet = union(el_merged$A, el_merged$B)
#'   return(gene_allnet)
#' }
#' 
get_allnodes = function(el_list){
  gene_allnet = sort(unique(unlist(el_list)))
  return(gene_allnet)
}