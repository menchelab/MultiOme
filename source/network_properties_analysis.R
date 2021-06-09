## Network Properties
## basic plots of nodes and edges
## Ize Buphamalai
## update Aug 25, 2020


# 0 - set up -------------
source("../functions/readdata_functions.R")
library(igraph)


# 1 - load required files -----------
g = process_graph_data("../data/network_edgelists/")

all_nodes = sapply(g, function(x) V(x)$name) %>% unlist %>% unique 

# citation count (# PMID per genes, queried and processed by INDRA)
citation_count = read_csv("../data/all_pmids_counts.csv", col_names = c("gene", "count"), col_types = "ci")
citation_count = citation_count %>% arrange(-count) %>% mutate(rank = 1:nrow(.))

# label of networks
network_details = read_tsv("../data/network_details.tsv")

citation_correlation = function(graph){
  #' take the citation count and correlate it with the degree of each graph
  
  degree_g = degree(graph)
  degree_g = tibble(gene = names(degree_g), degree = degree_g) 
  degree_g = inner_join(degree_g, citation_count, by = "gene") %>% arrange(rank) 
  #  degree_g =degree_g %>% filter(degree > 0) %>% group_by(gr=cut(degree, breaks=c(1,5,10,50,100,500,1000,5000), right=FALSE))
  
  #corval = cor(degree_g$count, degree_g$degree, method = "spearman")
  return(degree_g)
}

degree_citation_df = lapply(g, citation_correlation)

saveRDS(degree_citation_df, "../cache/degree_citation_correlation.RDS")

# 2 - compute properies -----------
print("computing topological properties")

# compute different network properties
n_nodes = sapply(g, vcount)
n_edges = sapply(g, ecount)
density = sapply(g, edge_density)
clustering = sapply(g, transitivity)
local_clustering = sapply(g, function(graph) transitivity(graph, type = "average"))
assortativity = sapply(g, function(graph) assortativity_degree(graph, directed = F))
social_bias = sapply(degree_citation_df, function(x) cor(x$count, x$degree, method = "spearman"))

# 3 - process the results and prepare for plotting

g_prop_df = tibble(network = names(n_nodes), n_nodes, n_edges, density, clustering, local_clustering, assortativity, social_bias)
g_prop_df = melt(g_prop_df, variable.name = "property")

# modifying the dataframe for plotting

g_prop_df=left_join(g_prop_df, network_details, by = "network") %>%
  mutate(
    # rename properties for readability
    property = factor(property, labels = c("Number of nodes", "Number of edges", "Edge density", "Global clustering", "Avg local clustering", "Assortativity", "Social bias")),
    # groups to plot
    group = ifelse(type == "Co-expression" & !is.na(type) , type, subtype),
    # change to factor
    group = factor(group, levels = rev(c("Co-essentiality",  "Biological Process", "Molecular Function", "Co-expression",  "Protein-protein interaction",
                                                     "Co-Pathway Membership", "Mammalian Phenotype", "Human Phenotype" ))),
    # alpha values (co-ex or not co-ex -- for aesthetic labelling)
    alphaval = ifelse(group == "Co-expression", "1", "0")
  )

saveRDS(g_prop_df, "../cache/network_complementarity_topological.RDS")






