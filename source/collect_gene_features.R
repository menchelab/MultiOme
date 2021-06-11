# Collecting gene features that can be used as predictor
# Ize Buphamalai
# Jun 2021


#' Features collected are
#' 1. Popularity - via citation counts
#' 2. # HPO terms matching genes
#' 3. # Associated GO terms
#' 4. # Associated reactome pathways


# Prediction based on popularity
citation_count = read_csv("../data/all_pmids_counts.csv", col_names = c("gene", "count"), col_types = "ci") %>% 
  arrange(-count) %>% 
  mutate(rank = 1:nrow(.)) %>%
  rename(PubMedRank = 'rank')

#patient_all_genes_with_features_df <- patient_all_genes_with_features_df %>%
#  left_join(., citation_count[,c('gene', 'PubMedRank')], by = c(Gene = 'gene'))

# prediction based on brain expression
gtex_avg_expression = read_tsv("../data/raw_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", skip = 2)

gtex_avg_expression_df <- tibble(Gene = gtex_avg_expression$Description, 
                                 AllBrainAvg = apply(gtex_avg_expression[,grepl("Brain", colnames(gtex_avg_expression))], 1, mean),
                                 FrontalCortex = gtex_avg_expression$`Brain - Frontal Cortex (BA9)`)

patient_all_genes_with_features_df <- patient_all_genes_with_features_df %>%
  left_join(., gtex_avg_expression_df)

# prediction based on pathway involvements
reactome_count_df <- read_tsv("../cache/reactome_pathway_counts.tsv", col_names = c("Gene", "n_pathways"), col_types = "ci", skip = 1)

patient_all_genes_with_features_df <- patient_all_genes_with_features_df %>%
  left_join(., reactome_count_df)

gene_feature_df <- full_join(citation_count, gtex_avg_expression_df, by = c(gene = "Gene")) %>%
  full_join(., reactome_count_df, by = c(gene = "Gene")) 

### add GO-based features
GO_count_by_genes <- read_tsv("../cache/GO_terms_count_per_gene.tsv")
gene_feature_df <- left_join(gene_feature_df, GO_count_by_genes)


### add HPO-based features
HPO_count_by_genes <- read_tsv("../cache/HPO_terms_count_per_gene.tsv")
gene_feature_df <- left_join(gene_feature_df, HPO_count_by_genes)

write_tsv(gene_feature_df, file = "../cache/gene_features.tsv")
