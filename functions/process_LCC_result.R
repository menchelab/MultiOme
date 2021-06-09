## function to process LCC results by computing pvals

# a function for processing LCC significance level for plotting
process_LCC_result = function(result_df, min_N_in_graph = 10, min_LCC_size = 5, network_annotate = T){
  #' @input result_df: object from LCC computation
  #' @input min_N_in_graph: minimum number of nodes in graph to include in the result
  #' @input min_LCC_size:  minimum number of LCC size to include in the result
  #' @output: processed_result_df: computed pvals 
  
  # process results
  processed_result_df = result_df %>% 
    mutate(LCC.pval = pnorm(LCC.zscore, lower.tail = F),
           network = factor(network),
           correctedPval = p.adjust(LCC.pval, method = "BH"),
           LCC.signif = cut(correctedPval, 
                            breaks=rev(c(1, 5e-2, 1e-2, 1e-3, 1e-4, 0)),
                            labels=rev(c("none", "*","**","***", "****")), ordered_result = T),
           minusLogpval = -log10(correctedPval)
    )
  
  # removed igraph object from the result
  if( "LCC.igraph" %in% colnames(processed_result_df)){
    processed_result_df = processed_result_df %>% select(-LCC.igraph) 
  }
  
  
  # filter for low quality lCC measures
  processed_result_df = processed_result_df %>% dplyr::filter(N_in_graph >= min_N_in_graph, LCC.size >= min_LCC_size)
  
  # annotate the network with network info
  if(network_annotate){
    #load network labels
    network_info = read_tsv("../data/network_details.tsv")
    
    # add network info
    processed_result_df = left_join(processed_result_df, network_info, by = "network")
    
    # count the amount of networks 
    sorted_networks_by_count = processed_result_df %>% ungroup %>% dplyr::filter(LCC.signif!="none")
    levels_networks = sorted_networks_by_count %>% count(subtype, sort = T) %>% pull(subtype) %>% as.character()
    levels_name = sorted_networks_by_count %>% count(name, sort = T)  %>% pull(name) %>% as.character()
    
    
    # reformat the result data frame with proper names
    
    processed_result_df = processed_result_df %>% 
      dplyr::filter(name %in% levels_name, subtype %in% levels_networks) %>%
      mutate(name = factor(name, levels = levels_name),
             subtype =  factor(subtype, levels = levels_networks))#, labels = network_abbr[levels_networks]))
    
  }
  
  return(processed_result_df)
}