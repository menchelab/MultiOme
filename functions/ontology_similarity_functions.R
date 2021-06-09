####
# set of functions related to ontology similarity

# Connect to server
#===============================

serverconnect = function(){
  require(dbplyr)
  require(DBI)
  db <- dbConnect(RMySQL::MySQL(),
                  dbname = "onto_based_gene_similarity",
                  user  = 'readonly',
                  password = 'ra4Roh7ohdee',
                  host = 'menchelabdb.int.cemm.at')
  return(db)
}

# Retrieve data from the server to store locally
#===============================
fetchDatafromServer = function(db, table, all_genes){
  similarity_data <- tbl(db, table)
  
  # take a subset of data (only for those that matters)
  # remove the ID columns
  sim_df = similarity_data %>% 
    filter(Name1 %in% all_genes & 
             Name2 %in% all_genes) %>% 
    select_at(vars(-contains("ID"))) %>% 
    collect()
  print(eval(sprintf("%s data fetched from the server", table)))
  
  return(sim_df)
}


# this function get all similarity scores between the dataframe and geneset provided
#===========================================
geneset_sim_measure_pair = function(sim_df, geneset1, geneset2){
  
  # the data might be at the first or the second column
  sim1  = sim_df %>% 
    filter(Name1 %in% geneset1 & 
             Name2 %in% geneset2) %>% 
    select_at(vars(-contains("Name")))
  
  sim2  = sim_df %>% 
    filter(Name1 %in% geneset2 & 
             Name2 %in% geneset1) %>% 
    select_at(vars(-contains("Name")))
  
  geneset_sim = rbind(sim1, sim2)
  mean_sim = colMeans(geneset_sim)
  
  return(mean_sim)
}


# this is the main calculation function that create table of similarity from selected tables
#===========================================
diseaseset_sim_measures = function(table_list, disgene_list){
  
  all_terms = names(disgene_list)
  all_genes = unique(unlist(disgene_list)) 
  
  # create an object for pairwise computation
  message("create combination of all terms")
  term_similarity_df = t(combn(all_terms, 2, simplify = T)) %>% as_tibble()
  colnames(term_similarity_df) = c("A","B")

  # objects for storing results
  disease_sim = list()
  random_sim = list()
  
  # connect to server
  db = serverconnect()
  
  # loop over different tables in the database
  for(table in table_list){
    # get the data
    #==============
    sim_df = fetchDatafromServer(db, table, all_genes)
    
    #  retrieve similarity scores
    #=============================
    message("Computing pairwise similarity...")
    disease_sim[[table]] = pbapply(term_similarity_df, 1, function(x)
      geneset_sim_measure_pair(sim_df, disgene_list[[x['A']]],
                               disgene_list[[x['B']]]))
    
    
    random_sim[[table]] = sim_df %>% select(-Name1, -Name2) %>% 
      summarise_all(funs(mean))
  }
  
  # bind terms to their similarity score
  disease_sim_df = do.call(rbind.data.frame, disease_sim)
  term_similarity_df = cbind(term_similarity_df, t(disease_sim_df))
  
  # add random similarity score
  random_sim_df = c(Name1="random", Name2 = "random", unlist(random_sim))
  term_similarity_df = rbind(term_similarity_df, random_sim_df)
  
  return(term_similarity_df)
}
