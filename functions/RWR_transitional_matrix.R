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


