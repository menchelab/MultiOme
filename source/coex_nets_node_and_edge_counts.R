# GTEx characterization for node and edge number at each filtering steps
# Ize Buphamalai
#
# prediction based on brain expression
gtex_avg_expression = read_tsv("../data/raw_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", skip = 2)

tissues_details <- read_csv("../data/tissue_combined_details.csv")

ensg_to_symbol <- readRDS("../cache/ensg_to_entrez_and_symbol_df.RDS")

protein_coding_sym <- read_tsv("../cache/coding_genes_in_gtex.tsv")
protein_coding_sym$protein_coding = T

ensg_to_symbol <- full_join(ensg_to_symbol, protein_coding_sym) %>%
  dplyr::rename("Description" = "symbol")

gtex_avg_expression_with_entrez <- full_join(ensg_to_symbol, gtex_avg_expression)


gtex_expression_count <- list()
gtex_expression_count$AllReported <- apply(gtex_avg_expression_with_entrez[,-c(1:8)],2, function(x) sum(x>0, na.rm = T))
gtex_expression_count$Protein_coding <- apply(gtex_avg_expression_with_entrez %>% filter(!is.na(protein_coding)) %>% select(-c(1:8))   ,2, function(x) sum(x>0, na.rm = T))
gtex_expression_count$avg_TPMover1 <- apply(gtex_avg_expression_with_entrez %>% filter(!is.na(protein_coding)) %>% select(-c(1:8))   ,2, function(x) sum(x>=1, na.rm = T))



source("../functions/readdata_functions.R")
library(igraph)
g = process_graph_data("../data/network_edgelists/")
coex_nets <- names(g)[grepl("coex", names(g))]
nodes <- lapply(coex_nets, function(x) V(g[[x]])$name)
names(nodes) <- coex_nets

gtex_expression_count$Final_nodes <- sapply(nodes, length)


# check with version before the disp filter
###########

count_nodes_from_edgelist <- function(el_file){
  el <- read_tsv(paste0("../cache/temp/new_coex_by_group//", el_file), col_types = "ff")
  #nodes <- length(union(levels(el$A), levels(el$B)))
  nodes <- length(union(levels(el$gene1), levels(el$gene2)))
  edges <- nrow(el)
  return(c(nodes = nodes, edges = edges))
}

files <- list.files("../cache/temp/new_coex_by_group//")

full_nodes = list()
for(i in files){
  print(i)
  full_nodes[[i]] = count_nodes_from_edgelist(i)
}

gtex_expression_count$BeforeDisparity <- sapply(full_nodes, function(x) x['nodes'])
gtex_expression_count$BeforeDisparity = ifelse(gtex_expression_count$BeforeDisparity < gtex_expression_count$Final_nodes, gtex_expression_count$Final_nodes, gtex_expression_count$BeforeDisparity)
###########

gtex_expression_count_df <-reshape2:: melt(gtex_expression_count)
colnames(gtex_expression_count_df) = c("count", "Measure")

gtex_expression_count_df <- gtex_expression_count_df %>%
  mutate(Measure = factor(Measure, levels =  c("AllReported",  "Protein_coding",  "BeforeDisparity", "Final_nodes", "avg_TPMover1"), label = c("All transcripts", "Protein\ncoding", "After\ncorr. cutoff", "After\ndisp. filter", "(Reference:\nTPM > 1)")))


edge_count <- list()
edge_count[["Before\ndisp. filter"]] <- sapply(full_nodes, function(x) x['edges'])
edge_count[["After\ndisp. filter"]] <- sapply(coex_nets, function(x) ecount(g[[x]]))
edge_count_df <-reshape2:: melt(edge_count)
colnames(edge_count_df) = c("count", "Measure")
edge_count_df$Measure <- factor(edge_count_df$Measure, levels = c("Before\ndisp. filter", "After\ndisp. filter"))

saveRDS(list(node_counts = gtex_expression_count_df, edge_counts = edge_count_df), "../cache/coex_nets_node_and_edge_counts.RDS")
