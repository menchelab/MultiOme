
# Process the localisation files and calculate z-score for all diseases and subnetworks
# for any gene sets and graph sets

###########
## load all requirements
############
library(pacman)
pacman::p_load("tidyverse", "igraph", "readr")
pacman::p_load_gh("psolymos/pbapply") # use apply functions with progress bars

source("~/Documents/projects/Multiome/functions/LCC_functions.R") 
source("~/Documents/projects/Multiome/functions/distance_functions.R") 
source("~/Documents/projects/Multiome/functions/readdata_functions.R") 

## load the arguments from the R arguments
##############

input = commandArgs(TRUE)
project_dir = input[1]
graph_path = input[2]   # input[1] -- expect a tsv file with 2 columns (edgelist) 
disgene_path = input[3] # input[2] -- expect a tsv file with 2 columns (disease, and ;-separated associated gene) 
suffix= input[4] # sub directory to store output
min_gene = ifelse(is.na(input[5]), 20, as.numeric(input[5])) # minimum amount of genes for a disease to consider in the measurement (default = 20) 
max_gene = ifelse(is.na(input[6]), 2000, as.numeric(input[6])) # maximum amount of genes for a disease to consider in the measurement (default = 2000) 
perform_disease_relationship = as.logical(input[7]) #logical, whether to perform disease-relationship part
distcal = as.logical(input[8])
## interactive settings
pbo = pboptions(type="txt")
#pboptions(pboptions(type = "timer"))

## check if input paths exists
for(paths in c(project_dir, graph_path, disgene_path)){
  if(!file.exists(paths)){
    name = deparse(paths)
    stop(sprintf("The path for %s: %s does not exist", deparse(paths), paths))
  }
}

# set working directory
setwd(project_dir)

# create cache and output directory
if(!file.exists("./cache")) system("mkdir cache")
if(!file.exists("./cache/distance_matrices")) system("mkdir cache/distance_matrices")
if(!file.exists("./output")) system("mkdir output")

# create the subdirectory for output
subdir = paste0("./output/" ,suffix)
if(!file.exists(subdir)) system(paste("mkdir", subdir))



## create list of graph 
###############
g = process_graph_data(graph_path)


## read disease gene-disease relationship file
##############
disgene = process_disease_genes_data(disgene_path, min_gene = min_gene, max_gene = max_gene)
disgene_pass_df = disgene$disgene_df
disgene_list = disgene$disgene_list

## Measuring LCC of each disease group in the network
####################################
# output file to store
out = sprintf("%s/LCC_and_distance_calculation_results.RDS", subdir)

#only perform if the output not exist
# 1. performing LCC calculation
if(!file.exists(out)){
  message('performing LCC calculation with random expectation')
  
  signif_df =  expand.grid(name = disgene_pass_df$name, network  = names(g), stringsAsFactors = F)
  
  LCC_df_cal = pbapply(signif_df, 1, function(x){
    LCC_randomisation_measure(graph = g[[x['network']]], 
                              nodesets = disgene_list[[x['name']]])})
  
  LCC_df = as_tibble(t(LCC_df_cal))
  
  # save igraph objects
  message('saving igraph objects of disease modules')
  
  LCC_igraph_cal = pbapply(signif_df, 1, function(x){
    LCC_igraph_compute(graph = g[[x['network']]], 
                       nodesets = disgene_list[[x['name']]])})
  
  
  signif_df = cbind(signif_df, LCC_df)
  signif_df$LCC.igraph = LCC_igraph_cal
  
  saveRDS(signif_df, out)
  
  message(paste('results saved as', out))
} else{
  message("LCC results already exist. Loading data")
  signif_df = readRDS(out)
}


# distance-based calculation
####################################
distance_exists = any(c("d_c", "d_s") %in% colnames(signif_df))

if(distcal & !distance_exists){
  message('performing distance calculation within geneset with random expectation')
  ## test for the signicance of random against random gene set
  
  # load distance matrix
  dismat = lapply(unique(signif_df$network), function(x) load_dismat(graph_name = x))
  
  dis_signif_results = pbapply(signif_df, 1, function(x){
    distance_cal(g = g[[x[['network']]]] ,dismat = dismat[[x[['network']]]], 
    #distance_cal(dismat = load_dismat(x[['network']]), 
                 vset1 = disgene_list[[x[['name']]]], silent = T)})
  
  names(dis_signif_results) = 1:length(dis_signif_results)
  
  dis_signif_results_df = t(as_tibble(dis_signif_results))  
  colnames(dis_signif_results_df) = names(dis_signif_results[[1]])
  
  
  # merge results and save data
  ####################################
  message('merging module LCC and distance results')
  signif_df = cbind(signif_df, dis_signif_results_df)
  saveRDS(signif_df, out)
}  


# determine whether disease relationship part (longer calculation) is required
if(perform_disease_relationship){
  message('performing disease relationship measurements')
  source("~/Documents/projects/Multiome/functions/main_disease_similarity_measure_function.R")
}

