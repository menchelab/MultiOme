#' Extracting HPO terms
#' Ize Buphamalai, CeMM
#' Update: June 2021
#' 

library(ontologyIndex)
data(hpo)

library(ontologySimilarity)

gene_hpo_terms_df <- read_tsv("../data/raw_data/HPO_phenotype_to_genes.txt", 
                              col_types = c('cc-c---'), col_names = c("HPOID", "HPOname", "gene"),
                              skip = 1) 

gene_HPO_terms <- split(gene_hpo_terms_df$HPOID, gene_hpo_terms_df$gene)

HPO_IC <- get_term_info_content(hpo, unique(gene_hpo_terms_df$HPOID)) 

# median
genes_ic5 <- sapply(gene_HPO_terms, function(x) sum(HPO_IC[x] >=  quantile(HPO_IC, 0.75), na.rm = T))

HPO_count_by_genes <- tibble(gene = names(gene_HPO_terms), 
                            HP_terms = sapply(gene_HPO_terms, length),
                            Informative_HP_terms = genes_ic5)


write_tsv(HPO_count_by_genes, "../cache/HPO_terms_count_per_gene.tsv")
