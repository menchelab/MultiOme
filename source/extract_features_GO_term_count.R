#' Extracting GO term counts per gene
#' Ize Buphamalai, CeMM
#' Update: June 2021
#' 

library(ontologyIndex)
data(go)

library(ontologySimilarity)
data(gene_GO_terms)
data(GO_IC)

#all_gene_names = citation_count$gene
#genes = gene_GO_terms[all_gene_names]
#genes = gene_GO_terms
#genes = genes[-which(sapply(genes, length)==0)]

bp <- go$id[go$name == "biological_process"]
genes_bp <- lapply(gene_GO_terms, function(x) intersection_with_descendants(go, roots=bp, x)) 
#data.frame(check.names=FALSE, `#terms`=sapply(genes, length), `#CC terms`=sapply(genes_bp, length))


mf <- go$id[go$name == "molecular_function"]
genes_mf <- lapply(gene_GO_terms, function(x) intersection_with_descendants(go, roots=mf, x)) 

#genes_bp_ic5 <- sapply(genes_bp, function(x) sum(GO_IC[x] > 5))

GO_count_by_genes <- tibble(gene = names(gene_GO_terms), 
                   BP_terms = sapply(genes_bp, length), 
                   MF_terms = sapply(genes_mf, length),
                   AllGO_terms = sapply(gene_GO_terms, length))


write_tsv(GO_count_by_genes, "../cache/GO_terms_count_per_gene.tsv")


# below is commented, for computing similarity matrix based on this
# define new similarity matrix by applying Resnik similarity instead of Lin, and use BMA for combining the similarity

#GO_similarity_matrix  <- get_sim_grid(ontology=go, information_content=GO_IC,
#                                      term_sets=genes_bp, term_sim_method = "resnik")
#write.csv2(GO_similarity_matrix, "./GO_similarity_mat.csv", quote = F)
#save(GO_similarity_matrix, file = "./GO_similarity_matrix_resnik.RData")