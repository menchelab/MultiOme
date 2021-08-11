# PPI subset properties analysis
# Ize Buphamalai
# Jul 2021

p_load(reshape2, patchwork, igraph)

source("../functions/readdata_functions.R")
g <- process_graph_data("../data/network_edgelist_additional//")
g$PPI_HIPPIE <- process_graph_data("../data/network_edgelists/ppi.tsv")[[1]]

# select  for plotting - those derived from PPI
g <- g[names(g)][grepl("PPI_HIPPIE|core", names(g))]


# code chunks from source/network_properties_analysis.R

all_nodes = sapply(g, function(x) V(x)$name) %>% unlist %>% unique 

# citation count (# PMID per genes, queried and processed by INDRA)
citation_count = read_csv("../data/all_pmids_counts.csv", col_names = c("gene", "count"), col_types = "ci")
citation_count = citation_count %>% arrange(-count) %>% mutate(rank = 1:nrow(.))

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

n_nodes = sapply(g, vcount)
n_edges = sapply(g, ecount)
density = sapply(g, edge_density)
clustering = sapply(g, transitivity)
local_clustering = sapply(g, function(graph) transitivity(graph, type = "average"))
assortativity = sapply(g, function(graph) assortativity_degree(graph, directed = F))
social_bias = sapply(degree_citation_df, function(x) cor(x$count, x$degree, method = "spearman"))

# 3 - process the results and prepare for plotting

g_prop_df = tibble(network = names(n_nodes), n_nodes, n_edges, density, clustering, local_clustering, assortativity, social_bias)
g_prop_df = reshape2::melt(g_prop_df, variable.name = "property")


g_prop_df =g_prop_df %>%
  mutate(
    # rename properties for readability
    property = factor(property, labels = c("Number of nodes", "Number of edges", "Edge density", "Global clustering", "Avg local clustering", "Assortativity", "Social bias")),
    # groups to plot
    type = ifelse(grepl("core", network), "Co-expression", "PPI"), 
    main_type = type,
    subtype = sapply(network, function(x) strsplit(x ,"_")[[1]][2]), 
    group = subtype, 
    # alpha values (co-ex or not co-ex -- for aesthetic labelling)
    alphaval = 0
  )


# get the stats from the original values

g_prop_df_orig <- readRDS( "../cache/network_complementarity_topological.RDS")
g_prop_coex <- g_prop_df_orig %>% filter(type == "Co-expression") %>% select(-source)
g_prop_coex$group = "Tissue-specific"

g_prop_Labelled_df <- rbind(g_prop_df, g_prop_coex)

# make HIPPIE goes up the rank
g_prop_Labelled_df$group <- fct_relevel(g_prop_Labelled_df$group, "HIPPIE", "HIPPIELargeScale", "HIPPIECurated", "Tissue-specific")
g_prop_Labelled_df$group <- fct_relevel(g_prop_Labelled_df$group, rev)
g_prop_Labelled_df$group <- factor(g_prop_Labelled_df$group,
                                   levels = c("core", "Tissue-specific", "HIPPIECurated", "HIPPIELargeScale", "HIPPIE"), 
                                   labels = c("Core (multi-tissue)", "Tissue-specific", "Curated", "Large scale", "Full PPI"))

saveRDS(g_prop_Labelled_df, "../cache/PPI_subset_and_CoexCore_properties.RDS")
