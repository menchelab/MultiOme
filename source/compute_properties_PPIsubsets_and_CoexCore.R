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

####
##  Compute citation count and local assortativity
#####
ppi_degree <- degree(g$PPI_HIPPIE)
ppi_degree = ppi_degree[ppi_degree>0]

p_k = degree_distribution(g$PPI_HIPPIE)
p_k = p_k[-1]

# q_k =  remain degree distribution
q_k = (1:max(ppi_degree) +1)*c(p_k[-1],0)/mean(ppi_degree)
#q_k = ppi_degree-1

mean_q  = sum(1:max(ppi_degree)*q_k)
variance_q = sum((mean(ppi_degree)-1:max(ppi_degree))^2 * q_k)


neighour_degree_avg = knn(g$PPI_HIPPIE)$knn 
neighour_degree_avg = neighour_degree_avg[names(ppi_degree)]

local_assortativity = (ppi_degree-1)*ppi_degree*(neighour_degree_avg - 1 - mean_q)/(2*ecount(g$PPI_HIPPIE)*variance_q)

# Load bioplex subsets
bioplex_baits_293T <- read_tsv("../data/raw_data/Bioplex_293T-baits.tsv")
bioplex_baits_293T$cell = "293T"
bioplex_baits_HCT116 <- read_tsv("../data/raw_data/Bioplex_HCT116-baits.tsv")
bioplex_baits_HCT116$cell <- "HCT116"

bioplex_baits <- rbind(bioplex_baits_293T, bioplex_baits_HCT116)

## compute degree and assortativity
cuts <- c(0,1,2,5, 10,20, 50,100, 200, 500,1000,5000)
citation_assortativity = full_join(citation_count,
                                   tibble(gene = names(ppi_degree), 
                                          local_assortativity = local_assortativity,
                                          degree = ppi_degree,
                                          neighbour_degree = neighour_degree_avg,
                                          degree_binned = cut(degree, breaks = cuts, labels = cuts[-1]),
                                          label = ifelse(local_assortativity < -0.002, gene, ""))) %>%
  filter(!is.na(local_assortativity)) %>%
  mutate(inBioPlexbait = gene %in% bioplex_baits$`Bait Symbol`)

# save cached data
saveRDS(list(citation_assortativity = citation_assortativity, 
             neighour_degree_avg = neighour_degree_avg, 
             ppi_degree = ppi_degree, 
             local_assortativity = local_assortativity) , "../cache/citation_assortativities_PPI.RDS")
