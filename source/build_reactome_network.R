## reactome graph building
## information retrived: 30 Jan 2019
## Ize Buphamalai
## Modified June 2021

library("igraph")
library(tidyverse)
library(parallel)
pacman::p_load(pbapply)
####  read pathway gene set from Reactome --------
reactome_gmt = scan("../data/raw_data/ReactomePathways.gmt", what = "", sep = "\n")
gmt_to_list = function(x){
  res = list()
  contents = strsplit(x, "\t")[[1]]
  res$name = contents[1]
  res$url = contents[2]
  res$genes = contents[3:length(contents)]
  return(res)
}


#### convert the association into list -----
reactome_pathways = lapply(reactome_gmt, gmt_to_list)

reactome = lapply(reactome_pathways, function(x) x[["genes"]])
names(reactome) = lapply(reactome_pathways, function(x) x[["name"]])

#remove pathways with only one gene (too specific)
reactome = reactome[!sapply(reactome, length)==1]
# remove pathways with over 300 genes (uninformative)
reactome = reactome[!sapply(reactome, length)>300]



######### build reactome network
d1 <- stack(reactome)

# take valid gene symbols
source("../functions/fn_source.R")
entrez_names <- IDconvert(unique(d1$values), from = "SYMBOL", to = "ENTREZID")

d1 <- d1 %>% filter(values %in% names(entrez_names)[!is.na(entrez_names)])

# create edge lists from d1
g_gene_pathway <- graph_from_data_frame(d1, directed = F)
adj_gene_pathway <- get.adjacency(g_gene_pathway, type = "both")

adj_gene_gene <- adj_gene_pathway %*% t(adj_gene_pathway)

gene_element <- colnames(adj_gene_gene) %in% unique(d1$values)

adj_gene_gene <- adj_gene_gene[gene_element, gene_element]

diag(adj_gene_gene) <- NA

gene_gene_df <-get.data.frame(graph_from_adjacency_matrix(adj_gene_gene, weighted = T, diag = F))

# store gene and pathway as id
gene_id <- unique(d1$values)
pathway_id <- unique(d1$ind)

# convert data into numeric
#d1$values <- factor(d1$values, levels = gene_id, labels = 1:length(gene_id) )
d1$ind <- factor(d1$ind, levels = pathway_id, labels = 1:length(pathway_id) )

reactome_by_element <- split(as.integer(d1$ind), d1$values)

shared_pathway <- pbapply::pbapply(gene_gene_df, 1, function(x) intersect(reactome_by_element[[x[1]]], reactome_by_element[[x[2]]]))#,
                                 #  cl = cl)

gene_gene_df$shared_pathways <- shared_pathway

# convert into df
reactome_member_df <- tibble(pathway = names(reactome), member = reactome) %>% 
  unnest(member) 

reactome_count_df <- reactome_member_df%>%
  count(member) %>%
  arrange(-n)


####  save reactome pathway counts -------
write_tsv(reactome_count_df, "../cache/reactome_pathway_counts.tsv")  
write_tsv(reactome_member_df, "../cache/reactome_pathway_association.tsv")  


# test_multiple scenarios of reactome graph build
g_reactome <- list()

for(i in 1:10){
  g_reactome[[i]] <- graph_from_data_frame(reactome_copathway_df[reactome_copathway_df$weight > i,], directed = FALSE)
}


{par(mfrow = c(5,2))

  for(i in 1:10){
    fit_power_law(g_reactome[[i]], mode = "ccdf")
  }
  
}

# select a cutoff of five

write_tsv(reactome_copathway_df[reactome_copathway_df$weight > 5,], "../data/network_edgelists/reactome_copathway.tsv")

#################
