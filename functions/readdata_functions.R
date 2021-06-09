## functions related to process graph and disease-gene association data


process_graph_data = function(graph_path, pattern = ".tsv", header = T, delim = "\t", col_types = "cc"){
  message("reading graph data")
  #if the input is a directory, read all the files in the directory
  if(grepl(pattern, graph_path)){
    graph_files = graph_path
  }
  else{
    graph_files = list.files(graph_path, pattern = pattern, full.names = T)
  }
  
  
  g = list()
  for(i in graph_files){
    # extract file name out of the directory path
    filename = strsplit(i, split = pattern, fixed = T)[[1]][1] #remove extension
    filename = strsplit(filename, "/")[[1]] # split all the slashes, and take the last one
    g_name = filename[length(filename)]
    
    #graph_el = read_tsv(i, comment = "#", col_types = "cc")
    graph_el = read_delim(i, comment = "#", col_names = header, delim = delim, col_types = col_types)
    
    #ensure that there is no NA elements
    graph_el = graph_el[!is.na(graph_el[,1]) & !is.na(graph_el[,2]),]
    
    graph = graph_from_data_frame(graph_el, directed = F)
    graph = igraph::simplify(graph) # remove multiples and loops
    g[[g_name]] = graph
    
    cmd = sprintf("read network %s: %i nodes, %i edges", g_name, vcount(graph), ecount(graph))
    message(eval(cmd))
  }
  return(g)
}


####################################

process_disease_genes_data = function(disgene_path, min_gene, max_gene){
  message("reading disease gene association data")
  
  disgene_df = read_tsv(disgene_path, comment = "#", col_types = "cc")
  colnames(disgene_df) = c("name", "gene")
  # split the ;-separated disease genes into a list of disease gene
  # and check for uniqueness
  disgene_df = disgene_df %>% rowwise() %>% mutate(genes_all = unique(strsplit(gene, split = ";")),
                                                   N = length(genes_all))
  
  # filter some diseases with too low or too high amount of genes. 
  # This was controlled by input parameters 4 and 5
  disgene_pass_df = disgene_df %>% dplyr::filter(N >= min_gene & N <= max_gene)
  
  # convert the disease gene association into list
  disgene_list = as.list.data.frame(disgene_pass_df$genes_all)
  names(disgene_list) = disgene_pass_df$name
  
  cmd = sprintf("read %i diseases, of total %i associated genes.", nrow(disgene_pass_df), length(unique(unlist(disgene_list))))
  print(eval(cmd))
  return(list(disgene_df = disgene_pass_df, disgene_list = disgene_list))
}
