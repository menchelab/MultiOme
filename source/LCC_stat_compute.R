# Precompute the mean and LCCs for all the networks
# Ize Buphamalai
# 20.04.2020

# This process will take ~30 minutes per network, if possible - use parallel computing 

source("./functions/readdata_functions.R")
source("./functions/LCC_functions.R")
source("./functions/LCC_randomisation_from_graph.R")
library(igraph)
g = process_graph_data("./data/network_edgelists/")

stat_compute_vals <- lapply(g, LCC_stat_compute)

saveRDS(stat_compute_vals, "./cache/LCC_stat_compute_vals_all_networks.RDS")
