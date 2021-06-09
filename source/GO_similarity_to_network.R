GosimBP <- readRDS("../cache/GO_similarity_matrix_biological_process.RDS")

#GosimBP[diag(GosimBP)] <- NA

#GosimBP[GosimBP == min(GosimBP, na.rm = T)] <- 0

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
