# make sure for each of this: el1 is always less than el2, and no NAs
# Ize Buphamalai, last updated 17.07.19
process_edgelist = function(el, distinct = T){
  #' @input: a list of named edgelist
  #' @output: processed edgelist

  # remove edges with at least one node without label
  el = el[!is.na(el[,1]) & !is.na(el[,2]),]
  
  # remove edges where two nodes are the same
  el = el[el[,1] != el[,2],]

  swaprows = as.vector(el[,1] > el[,2])
  el_sorted = el
 # el_sorted[,1][swaprows] = el[,2][swaprows]
#  el_sorted[,2][swaprows] = el[,1][swaprows]
  el_sorted[swaprows,1] = el[swaprows,2]
  el_sorted[swaprows,2] = el[swaprows,1]
  if(distinct){
    el_sorted = el_sorted %>% distinct_at(vars(1,2))
  }
  return(el_sorted)
}