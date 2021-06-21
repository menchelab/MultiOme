GosimBP <- readRDS("../cache/GO_similarity_matrix_biological_process.RDS")

# similarity among rows

GOsimBP_mean <- apply(GosimBP, 1, function(x) mean(x, na.rm = T))
GOsimBP_median <- apply(GosimBP, 1,  function(x) median(x, na.rm = T))

# Load functions needed
source("../functions/sim_mat_functions.R")

#Remove the minimally relevant similarity information
GOBP_sim_mat_filter = remove_minval(GosimBP, n=1)


GOsimBP_mean <- apply(GOBP_sim_mat_filter, 1, function(x) mean(x, na.rm = T))
GOsimBP_median <- apply(GOBP_sim_mat_filter, 1,  function(x) median(x, na.rm = T))


#hist(k, breaks = 100, main = "Non-weighted degree distribution before filtering", col = "lightblue")
g_unfiltered <- graph_from_adjacency_matrix(GOBP_sim_mat_filter>5, mode = "undirected")

{par(mfrow = c(2,1))
  fit_power_law(g_unfiltered, mode = "ccdf")
  fit_power_law(g_unfiltered, mode = "pdf")}

paste0("nodes: ",length(V(g_unfiltered)))
paste0("edges: ",length(E(g_unfiltered)))


# now apply the backbone filtering algorithm
GO_pvalmat = backbone_cal(GOBP_sim_mat_filter)  


#####################################################################
# Apply pval = 0.05 to make adjacency matrix
#####################################################################

adjacency_GO = adj_mat_from_pvalmat(GO_pvalmat, threshold = 0.05)

g_GO_p05 = graph_from_adjacency_matrix(adjacency_GO, mode = "undirected")

paste0("nodes: ",length(V(g_GO_p05)))
paste0("edges: ",length(E(g_GO_p05)))


{par(mfrow = c(2,1))
  fit_power_law(g_GO_p05, mode = "ccdf")
  fit_power_law(g_GO_p05, mode = "pdf")}

######
GOBP_sim_mat_filter_selected <- GOBP_sim_mat_filter[rownames(adjacency_GO), colnames(adjacency_GO)]
GOBP_sim_mat_filter_dispfilter <- GOBP_sim_mat_filter_selected*adjacency_GO
#GOBP_sim_mat_filter_dispfilter[GOBP_sim_mat_filter_dispfilter==NA] <- 0


GOsimBP_dispfilter_median <- apply(GOBP_sim_mat_filter_dispfilter, 1,  function(x) median(x, na.rm = T))


# plot similarity histogram before vs after dispfilter ----------------
sim_score_file <- "../cache/GO_similarity_score_binned_before_and_after_dispfilter.tsv"
if(!file.exists(sim_score_file)){
  all_int_cut_df <- cut(GOBP_sim_mat_filter , breaks = seq(0.5,12,0.5), labels = 1:23)
  filtered_int_cut_df <- as.numeric(all_int_cut_df) * as.vector(GO_pvalmat < 0.05)
  
  all_int_cut <- table(all_int_cut_df)
  filtered_int_cut <- table(filtered_int_cut_df)
  
  sim_score_df <- tibble(group = 1:23, min_score = seq(0.5,11.5, 0.5))
  sim_score_df <- sim_score_df %>%
    mutate(all_int = all_int_cut[as.character(group)],
           remain_int = filtered_int_cut[as.character(group)])
  
  sim_score_df$all_int[is.na(sim_score_df$all_int)] = 0
  sim_score_df$remain_int[is.na(sim_score_df$remain_int)] = 0
  
  sim_score_plot_df <- sim_score_df %>%
    pivot_longer(cols = c(all_int, remain_int), names_to = "type", values_to = "count") %>%
    mutate(count = as.numeric(count)/2) %>%
    filter(count > 0)
  
  write_tsv(sim_score_plot_df, file = sim_score_file)
} else{
  sim_score_plot_df <- read_tsv(sim_score_file)
}

p <- ggplot(sim_score_plot_df, aes(x = min_score, y = count, fill = type)) + 
  geom_col(position = position_dodge2(preserve = "single")) + 
  scale_y_log10() +
  theme_cowplot() +
  ylab("# Interactions") +
  xlab("Similarity score") +
  scale_fill_manual(values = c("grey70", "salmon"))

ggsave("../Figs/barplot_GO_network_before_and_after_dispfilter.pdf", width = 6, height = 2.5)
  
