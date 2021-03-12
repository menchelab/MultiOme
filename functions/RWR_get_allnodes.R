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