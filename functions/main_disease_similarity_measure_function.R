#' compute topological separation between modules
#' Note: this code need to run after the main_localisation code 
###################################

filter = dplyr::filter

# DISTANCE AMONG DISEASE GENE PAIRS
###################################
# create unique combination of diseases
DisDis_df = t(combn(disgene_pass_df$name, 2, simplify = T)) %>% as_tibble()
colnames(DisDis_df) = c("Name.x", "Name.y")
# create custom index
DisDis_df$index = 1:nrow(DisDis_df)

# create full combination of disease pairs and network, using index to merge
DisDisNetwork_df = expand.grid(index = 1:nrow(DisDis_df), network = names(g))

# merge the data
DisDisNetwork_df = full_join(DisDisNetwork_df, DisDis_df, by = "index") %>%
  # get the N_in_graph for disease A
  left_join(., signif_df[,c("name", "network", "N_in_graph", "d_c_mean", "d_s_mean")], by = c("network", "Name.x" = "name")) %>%
  # get the N_in_graph for disease B
  left_join(., signif_df[,c("name", "network", "N_in_graph", "d_c_mean", "d_s_mean")], by = c("network", "Name.y" = "name"))

DisDisNetwork_df_filtered = DisDisNetwork_df %>% filter(N_in_graph.x >= 10 & N_in_graph.y >= 10)

DisDisNetwork_results = pbapply(DisDisNetwork_df_filtered, 1, function(x){
  distance_cal(dismat = dismat[[x[['network']]]], 
  #distance_cal(dismat = load_dismat(x[['network']]), 
               vset1 = disgene_list[[x[['Name.x']]]],
               vset2 = disgene_list[[x[['Name.y']]]],silent = T)})

names(DisDisNetwork_results) = 1:length(DisDisNetwork_results)

DisDisNetwork_results_df = t(as_tibble(DisDisNetwork_results))  
colnames(DisDisNetwork_results_df) = names(DisDisNetwork_results[[length(DisDisNetwork_results)]])

# merge results and save data
####################################
message('merging disease pair comparison results')
DisDisNetwork_df_filtered = cbind(DisDisNetwork_df_filtered,  DisDisNetwork_results_df)

# Compute the separation
DisDisNetwork_df_filtered = DisDisNetwork_df_filtered %>% as_tibble %>%
  mutate(Separation_s = unlist(d_s_mean) - (unlist(d_s_mean.x) + unlist(d_s_mean.y))/2,
         Separation_c = unlist(d_c_mean) - (unlist(d_c_mean.x) + unlist(d_c_mean.y))/2
  )

out = paste0(subdir,"/pairwise_distance_calculation_results.RDS")
saveRDS(DisDisNetwork_df_filtered, out)

message(paste('results saved as', out))