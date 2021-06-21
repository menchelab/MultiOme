# VIsualising GO tree structure 

library(ontologyIndex)
data(go)


library(ontologySimilarity)
data(gene_GO_terms)
data(GO_IC)


GO_feature_df <-  tibble(name = go$id, description = go$name, ic = GO_IC[go$id]) %>%
  filter(!name %in% go$obsolete, grepl("GO:", name), !is.na(ic))


go_df <- stack(go$children)
colnames(go_df) <-c("to","from")

go_df <- go_df[,c(2,1)]
go_df <- go_df %>% 
  mutate(from = as.character(from),
         to = as.character(to)) %>%
  filter(from %in% GO_feature_df$name, to %in% GO_feature_df$name) %>%
  mutate(ic_from = GO_IC[from], 
         ic_to = GO_IC[to],
         ic_diff = ic_to - ic_from) %>%
  distinct(to, .keep_all = T)


### Load GO data for detailed inspection

# CO sim mat
#GosimBP <- readRDS("../cache/GO_similarity_matrix_biological_process.RDS")

#{sim_allgenes <- GosimBP[lower.tri(GosimBP)]

#{pdf("../Figs/histogram_GOBP_similarity.pdf", height = 3.35, width = 4.75)
#  hist(sim_allgenes, breaks = 100, col = "#F8B100", main = "GOBP pairwise similarity", xlab = "Similarity score", border = "white")
#  dev.off()
#}}


# Load functions needed
source("../functions/sim_mat_functions.R")

#Remove the minimally relevant similarity information
GOBP_sim_mat_filter = remove_minval(GosimBP, n=1)

# now apply the backbone filtering algorithm
GO_pvalmat = backbone_cal(GOBP_sim_mat_filter)  

# =======

# functions creating sets of plots for detailed inspection 

GO_pairwise_compare_plot = function(gene1, gene2){
  require(ggraph)
  require(tidygraph)
  
  terms1 <- gene_GO_terms[[gene1]]
  terms2 <- gene_GO_terms[[gene2]]
  
  terms <- union(terms1, terms2)
  
  all_ancestors <- unlist(lapply(terms, function(x) go$ancestors[[x]])) 
  
  
  all_ancestors <- get_ancestors(go, terms)
  
  bp_ancestors <- setdiff(intersect(all_ancestors, get_descendants(go, roots = "GO:0008150")), "GO:0008150")
  
  go_dat <- list()
  
  go_dat$edges <- go_df %>% filter(from %in% bp_ancestors | to %in% bp_ancestors) #%>% 
  # filter(!(from %in% terms & !to %in% terms))# | to %in% bp_ancestors)
  
  
  go_dat$vertices <- GO_feature_df %>% 
    filter(name %in% unique(c(go_dat$edges$from, go_dat$edges$to))) %>% 
    mutate(interm1 = name %in% terms1,
           interm2 = name %in% terms2,
           label = ifelse(interm1 & interm2, "both", ifelse(interm1, gene1, ifelse(interm2, gene2, "none"))),
           label = factor(label, levels = c("both", gene1, gene2, "none"))
    ) 
  
  graph <- tbl_graph(go_dat$vertices, go_dat$edges)
  
  
  # this one is good
  #ggraph(graph, 'tree') + 
  #  geom_edge_diagonal(edge_colour = "grey") +
  #  geom_node_point(col = "red", aes(alpha =interm1))+
  #  scale_alpha_manual(values = c(0,1))
  
  p_go_dendro <- ggraph(graph, 'tree') + 
    geom_edge_diagonal(edge_colour = "grey") +
    geom_node_point(aes(alpha = interm1 | interm2, col = label))+
    #geom_node_text(aes( label = ifelse(interm1 & interm2, description, ""))) +
    geom_node_text(aes( label = ifelse(ic > 9.85 & interm2, description, ""))) +
    scale_alpha_manual(values = c(0,1)) +
    scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb',"#FFFFFF")) +
    guides(alpha = F) +
    theme_nothing()
  
  ggsave(sprintf("../Figs/GO_tree_label_%s_%s.pdf", gene1, gene2), p_go_dendro, width = 9, height = 3)  
  
  # # this one is also good 
  #ggraph(graph, 'dendrogram', circular = TRUE) + 
  #  geom_edge_elbow(edge_colour = "grey") +
  #  geom_node_point(aes(alpha = interm1 | interm2, col = label))+
  #  scale_alpha_manual(values = c(0,1)) +
  #  scale_color_manual(values = c('#66c2a5','#fc8d62','#8da0cb',"#FFFFFF")) +
  #  guides(alpha = F) +
  #  theme_nothing() +
  #  coord_fixed()
  
  #information_content <- descendants_IC(go)
  
  
#  terms1 <- go_dat$vertices %>% arrange(label) %>% filter(interm1) %>% pull(name)
#  terms2 <- go_dat$vertices %>% arrange(label) %>% filter(interm2)  %>% pull(name)
  
#  sim_mat <- ontologySimilarity::get_term_sim_mat(ontology=go,method = "resnik", row_terms = terms1, col_terms = terms2, information_content = GO_IC)
  
 # heatmap(sim_mat, Rowv = NA, Colv = NA, col = RColorBrewer::brewer.pal(11,'RdYlBu'))     
  
  
  common_terms_df <- go_dat$vertices %>% filter(label != "none") %>% select(-interm1, -interm2)
  
  p_bar <- ggplot(common_terms_df %>% count(label)) + geom_bar(aes(x = "count", y = n, fill = label), stat="identity") +
    scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb')) +
    ylab("# GO terms") +
    xlab("") +
    guides(fill = F) +
    theme_cowplot()
  ggsave(sprintf("../Figs/GO_stat_bar_%s_%s.pdf", gene1, gene2), p_bar, width = 1, height = 3)  
  
  
  p_density <- ggplot(common_terms_df) + geom_histogram(aes(x = ic, fill = label)) +
    facet_grid(label ~ .) +
    scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb')) +
    ylab("Terms") +
    xlab("Information Content") +
    guides(fill = F) +
    theme_cowplot()
  
  p_GO_stat <- p_bar + p_density
  ggsave(sprintf("../Figs/GO_stat_%s_%s.pdf", gene1, gene2), p_density, width = 3, height = 3)  
  
  
  sim_df <- tibble(gene = rownames(GO_pvalmat),
                   sim_gene1 = GOBP_sim_mat_filter[gene1,], 
                   sim_gene2 = GOBP_sim_mat_filter[gene2,],
                   pval_gene1 =GO_pvalmat[gene1,],
                   pval_gene2 =GO_pvalmat[gene2,]) %>%
    mutate(overthreshold1 = -log10(pval_gene1) >= -log10(0.05),
           overthreshold2 = -log10(pval_gene2) >= -log10(0.05))
  
  
  
  p_scatter_gene1 <- ggplot(sim_df %>% filter(!is.na(overthreshold1)),  aes(x = sim_gene1, 
                                                                            y = -log10(pval_gene1), 
                                                                            label = ifelse(overthreshold1, gene, ""))) + 
    geom_point(aes(col = overthreshold1, alpha = overthreshold1)) +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    ggrepel::geom_label_repel() + theme_cowplot() +
    scale_color_manual(values = c("grey80", "grey20")) +
    scale_alpha_manual(values = c(0.25,0.75)) +
    guides(color = F, alpha = F) +
    xlab("Similarity score") +
    ylab("Disparity p-value") +
    ggtitle(sprintf("Similarity score: %s", gene1))
  
  p_scatter_gene2 <- ggplot(sim_df %>% filter(!is.na(overthreshold2)),  aes(x = sim_gene2, 
                                                                            y = -log10(pval_gene2), 
                                                                            label = ifelse(overthreshold2, gene, ""))) + 
    geom_point(aes(col = overthreshold2, alpha = overthreshold2)) +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    ggrepel::geom_label_repel() + theme_cowplot() +
    scale_color_manual(values = c("grey80", "grey20")) +
    scale_alpha_manual(values = c(0.25,0.75)) +
    guides(color = F, alpha = F) +
    xlab("Similarity score") +
    ylab("Disparity p-value") +
    ggtitle(sprintf("Similarity score: %s", gene2))
  
  
  p_scatter_selected <- p_scatter_gene1 + p_scatter_gene2
  ggsave(sprintf("../Figs/scatter_plot_similarity_score_%s_%s.pdf", gene1, gene2), p_scatter_selected, width = 9, height = 4)
}

# explore examples of different gene pairs

gene_pairs_to_plot <- tibble(gene1 = c("TP53", "SLC12A9","OR10J3","EMC3"),
                             gene2 = c("TGFB1", "ABCC3", "OR13G1", "EMC1"))


# plot and save results for these gene pairs
apply(gene_pairs_to_plot, 1, GO_pairwise_compare_plot)

# get data frame for shared genes


gene_pairs_to_plot$term_union <- apply(gene_pairs_to_plot, 1, function(x) length(unique(unlist(gene_GO_terms[x]))))
gene_pairs_to_plot$term_gene1 <- sapply(gene_GO_terms[gene_pairs_to_plot$gene1], length)
gene_pairs_to_plot$term_gene2 <- sapply(gene_GO_terms[gene_pairs_to_plot$gene2], length)
gene_pairs_to_plot <- gene_pairs_to_plot %>%
  mutate(term_overlap = term_gene1 + term_gene2 - term_union,
         term_gene1_only = term_gene1 - term_overlap,
         term_gene2_only = term_gene2 - term_overlap)

gene_pairs_to_plot_df <- gene_pairs_to_plot %>%
  select(gene1, term_gene1_only, term_overlap, term_gene2_only) %>%
  pivot_longer(cols = c(term_gene1_only, term_overlap, term_gene2_only), names_to = "set", values_to = "n") %>%
  mutate(set = factor(set, levels = c("term_gene1_only", "term_overlap", "term_gene2_only"))) %>%
  group_by(gene1) %>%
  mutate(frac_n = n/sum(n))
         
p_bar_frac <-ggplot(gene_pairs_to_plot_df) + geom_bar(aes(x = gene1, y = frac_n, fill = set), stat="identity") +
  scale_fill_manual(values = c('#fc8d62','#66c2a5','#8da0cb')) +
  ylab("#fraction GO terms") +
  xlab("") +
  guides(fill = F) +
  theme(axis.text.x = element_blank()) +
  theme_cowplot()


p_bar_union  = ggplot(gene_pairs_to_plot) + geom_col(aes(x = gene1, y = term_union), fill = "grey50") +
  ylab("#GO terms") +
  xlab("") +
  scale_y_log10() +
  guides(fill = F) +
  theme(axis.text.x = element_blank()) +
  theme_cowplot()


ggsave(sprintf("../Figs/GO_stat_bar_all_example.pdf"), p_bar_frac, width = 2, height = 2.5)  
ggsave(sprintf("../Figs/GO_stat_bar_union_all_example.pdf"), p_bar_union, width = 2, height = 1.5)  

# Global similarity score
sim_df <- tibble(sim = GOBP_sim_mat_filter[lower.tri(GOBP_sim_mat_filter)],
                 pval =GO_pvalmat[lower.tri(GO_pvalmat)])

sim_df_plot <- sim_df %>% filter(sim > 5)


p <- ggplot(sim_df_plot) + geom_hex(aes(sim, pval), bins = 100) + 
  scale_fill_gradient(low = "lightblue", high = "darkblue") + 
  xlab("Similarity score") +
  ylab("Disparity p-value") +
  scale_y_log10() +
  theme_cowplot()
ggsave("../Figs/GO_sim_and_pval.pdf", p, width = 5, height = 3)

p <- ggplot(term_df) + 
  geom_point(aes(BP_terms, degree, col = log10(count))) + 
  scale_color_viridis_c() + scale_x_log10() + scale_y_log10() +
  theme_cowplot()+
  xlab("# GO terms (BP)") +
  ylab("GOBP degree")
ggsave("../Figs/GO_degree_and_term_count.pdf", p, width = 4.5, height = 3.75)



