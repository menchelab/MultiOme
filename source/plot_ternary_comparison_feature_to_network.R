# Make ternary plot for three different angles in turning features into network
# social bias, its feature and its degree

# perform the test in BP, HP, reactome and brain expression


# load required networks
g <- process_graph_data(paste0("../data/network_edgelists/", c("coex_BRO", "GOBP", "GOMF", "HP", "reactome_copathway"), ".tsv"))

degrees <- lapply(g, function(x) tibble(degree = degree(x), gene = names(degree(x))))
names(degrees) <- paste0("degree_", names(degrees))
degree_df <- bind_rows(degrees, .id = "network") %>%
  pivot_wider(., names_from = network, values_from = degree)

#
gene_feature_df <- read_tsv(file = "../cache/gene_features.tsv")

term_df <- left_join(gene_feature_df, degree_df) 



pacman::p_load(Ternary)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf("../Figs/ternary_plot_curated_databases.pdf", width = 12, height = 3.5)
{par(mfrow = c(1, 4), mar = rep(0.2, 4))
  
  # GOBP --------
  coordinates <- term_df %>%select(gene, count, degree_GOBP, BP_terms) %>% arrange(-BP_terms)
  
  coordinates <- coordinates[!apply(coordinates[,2:4],1, function(x) any(is.na(x))),]
  top_genes <- coordinates$gene[1:10]
  names(data_points) <- coordinates$gene[1:10]
  coordinates <- apply(-coordinates[,2:4], 2, rank)
  
  #data_points <- lapply(1:length(top_genes), function(x) as.double(coordinates[x,]))
  #names(data_points) <- top_genes
  
  TernaryPlot(#atip = 'PubMed count', btip = 'Degree', ctip = 'GO terms',
    #alab = "\u2190 More PubMed count", blab = "\u2190 more network degree", clab = "More GO terms associated  \u2192",
    axis.labels = seq(1, nrow(coordinates), by = round(nrow(coordinates)/10)))
  ColourTernary(TernaryDensity(coordinates, resolution = 10L))
  TernaryPoints(coordinates, col = add.alpha("red", alpha=0.1), pch = '.')
  TernaryDensityContour(coordinates, resolution = 30L)
  #data_points <- list(TP53 = coordinates[1,])
  #AddToTernary(text, data_points, names(data_points), cex = 0.8, font = 2)
  
  
  # HPO ---------
  coordinates <- term_df %>%select(count, degree_HP, HP_terms)
  
  coordinates <- coordinates[!apply(coordinates,1, function(x) any(is.na(x))),]
  
  coordinates <- apply(-coordinates, 2, rank)
  
  TernaryPlot(#alab = "\u2190 More PubMed count", blab = "\u2190 more network degree", clab = "More phenotypes associated  \u2192",
    axis.labels = seq(1, nrow(coordinates), by = round(nrow(coordinates)/10)))
  ColourTernary(TernaryDensity(coordinates, resolution = 10L))
  TernaryPoints(coordinates, col = add.alpha("red", alpha=0.1), pch = '.')
  TernaryDensityContour(coordinates, resolution = 30L)
  # data_points <- list(TP53 = coordinates[1,])
  # AddToTernary(text, data_points, names(data_points), cex = 0.8, font = 2)
  
  
  # reactome ---------
  coordinates <- term_df %>%select(count, degree_reactome_copathway, n_pathways)
  
  coordinates <- coordinates[!apply(coordinates,1, function(x) any(is.na(x))),]
  
  coordinates <- apply(-coordinates, 2, rank)
  
  TernaryPlot(#alab = "\u2190 More PubMed count", blab = "\u2190 more network degree", clab = "More pathways associated  \u2192",
    axis.labels = seq(1, nrow(coordinates), by = round(nrow(coordinates)/10)))
  ColourTernary(TernaryDensity(coordinates, resolution = 10L))
  TernaryPoints(coordinates, col = add.alpha("red", alpha=0.1), pch = '.')
  TernaryDensityContour(coordinates, resolution = 30L)
  #data_points <- list(TP53 = coordinates[1,])
  #AddToTernary(text, data_points, names(data_points), cex = 0.8, font = 2)
  
  
  # Brain -----------
  coordinates <- term_df %>%select(count, degree_coex_BRO, AllBrainAvg)
  
  coordinates <- coordinates[!apply(coordinates,1, function(x) any(is.na(x))),]
  
  coordinates <- apply(-coordinates, 2, rank)
  
  TernaryPlot(#alab = "\u2190 More PubMed count", blab = "\u2190 more network degree", clab = "higher expression level  \u2192",
    axis.labels = seq(1, nrow(coordinates), by = round(nrow(coordinates)/10)))
  ColourTernary(TernaryDensity(coordinates, resolution = 10L))
  TernaryPoints(coordinates, col = add.alpha("red", alpha=0.1), pch = '.')
  TernaryDensityContour(coordinates, resolution = 30L)
  #  data_points <- list(TP53 = coordinates[1,])
  #  AddToTernary(text, data_points, names(data_points), cex = 0.8, font = 2)
  
}
dev.off()


