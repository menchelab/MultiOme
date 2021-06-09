## reactome graph building
## information retrived: 30 Jan 2019
## Ize Buphamalai
## Modified June 2021


library("igraph")

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

# convert into df
reactome_member_df <- tibble(pathway = names(reactome), member = reactome) %>% 
  unnest(member) 

reactome_count_df <- reactome_member_df%>%
  count(member) %>%
  arrange(-n)


####  save reactome pathway counts -------
write_tsv(reactome_count_df, "../cache/reactome_pathway_counts.tsv")  
write_tsv(reactome_member_df, "../cache/reactome_pathway_association.tsv")  



#### converting into networks --------

# make pair combn of genes in the same pathways
reactome_list_df = lapply(reactome, function(x) as_tibble(t(combn(x,2, simplify = T))))
names(reactome_list_df) = names(reactome)

# combine them into one
reactome_df = bind_rows(reactome_list_df)


# count each instance
reactome_copathway_df = reactome_df %>% 
  group_by(V1, V2) 

  colnames(reactome_copathway_df) = c("gene1", "gene2", "weight")

## remove the pair with less than four pathways
#reactome_copathway_filtered_df = reactome_copathway_df[reactome_copathway_df$weight > 1,]
reactome_copathway_filtered_df = reactome_copathway_df


# apply disparity filter (optional scenario)
#source("~/Documents/projects/universal_fn/network_functions.R")

#weighted adj mat from all data
#g_reactome_full = graph_from_data_frame(reactome_copathway_filtered_df, directed = F)

#g_reactome_full = graph_from_data_frame(reactome_copathway_filtered_df, directed = F)
#adj_reactome_full = as_adjacency_matrix(g_reactome_full, type = "both", attr = "weight")
#disp_filter = backbone_cal(adj_reactome_full)

#adj_filter = disp_filter < 0.05
#g_reactome_disp_filter = graph_from_adjacency_matrix(adj_filter, mode = "undirected")
#el_greactome_disp_filter = as_edgelist(g_reactome_disp_filter) %>% as_tibble() 

#write_tsv(el_greactome_disp_filter, "./network_edgelists/reactome_copathway_disp_filter.tsv")



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
