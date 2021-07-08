## Preparing additional PPI networks for analyses
## Ize Buphamalai

source("../functions/fn_source.R")
source("../functions/readdata_functions.R")


library(pacman)

# HuRI- process networks using gene names as identifiers
HuRI_df <- read_tsv("../data/raw_data/HuRI.tsv", col_names = c("A","B"))
HuRI_df$A <- IDconvert(HuRI_df$A, from = "ENSEMBL", to = "SYMBOL")
HuRI_df$B <- IDconvert(HuRI_df$B, from = "ENSEMBL", to = "SYMBOL")

HuRI_df <- distinct(HuRI_df,A,B)
write_tsv(HuRI_df, "../data/network_edgelist_additional/ppi_HuRI.tsv")

# BioPlex
Bioplex_df <- read_tsv("../data/raw_data/BioPlex.tsv")
write_tsv(Bioplex_df %>% select(SymbolA, SymbolB), "../data/network_edgelist_additional/ppi_BioPlex.tsv")

# core transcriptional modules
coex_el_sum <- readRDS("../cache/coexpression_raw_edge_counts.RDS")
coex_core <- coex_el_sum %>% filter(n>5) %>% select(A,B)
write_tsv(coex_core, file = "../data/network_edgelist_additional/coex_core.tsv")



p_load(igraph)

g <- list()

# Add PPI
g$PPI_HIPPIE <- process_graph_data("../data/network_edgelists/ppi.tsv")[[1]]

g$ppi_HuRI <- process_graph_data("../data/network_edgelist_additional/ppi_HuRI.tsv")[[1]]
g$ppi_BioPlex <- process_graph_data("../data/network_edgelist_additional/ppi_BioPlex.tsv")[[1]]

# measure overlaps between these data
g$PPI_LargeScale <-  g$ppi_HuRI %u% g$ppi_BioPlex
write_tsv(as_data_frame(g$PPI_LargeScale), file = "../data/network_edgelist_additional/PPI_LargeScale.tsv")

g$PPI_HIPPIECurated <- g$PPI_HIPPIE %m% g$PPI_LargeScale
write_tsv(as_data_frame(g$PPI_HIPPIECurated), file = "../data/network_edgelist_additional/PPI_HIPPIECurated.tsv")

# interrsection
g$PPI_HIPPIELargeScale <- g$PPI_HIPPIE %s% g$PPI_LargeScale
write_tsv(as_data_frame(g$PPI_HIPPIELargeScale), file = "../data/network_edgelist_additional/PPI_HIPPIELargeScale.tsv")