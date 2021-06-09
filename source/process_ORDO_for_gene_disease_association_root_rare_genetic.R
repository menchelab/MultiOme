# Process the gene association data from Orphanet and map them to the top branches on the ORDO

# Pisanu Ize Buphamalai
# Update: 9 April 2020

# Re analyse Orphanet data: making 'Rare genetic disease' becoming a root (ORPHA:98053)
# Observation: the number of genes for some disease groups were substantially decreased


###############
# load required packages and data
###############
#Read orphanet gene association data (immediate association)
source("read_orphanet_gene_association_data.R")

# Load the orphanet gene disease association (with ancestor terms)
ExtendedOrphaAssociationFile = "../cache/orphanet_disease_gene_association_with_ancestors.RDS"
if(!file.exists(ExtendedOrphaAssociationFile)){
  source("../source/Orphanet_annotate_genes_to_ancestors.R")
} else{
 combinedDiseaseGeneList =  readRDS(ExtendedOrphaAssociationFile)
}

###############
# create a function to retrieve associated genes given a root term
###############

retrieve_all_gene_association = function(root_term, return_genelist = T){
  # if not return_genelist, disease descendant list will be return
  if(!exists("ordo")){
    library(pacman)
    pacman::p_load("rols", "stringr")
    
    ordo = Ontology("ordo")
  }
  # get the ordo object from the term
  root_disease = term(ordo, root_term)
  
  # retrieve all children of the roots (the disease groups to collect associated genes from)
  terms_disease = children(root_disease)
  
  # get the IDs of disease groups
  terms_disease_id = termId(terms_disease)
  
  # compute descendants of all groups
  #disease_descendants = lapply(terms_disease, descendants)
  disease_descendants <- lapply(terms_disease_id, function(x) descendants(term(ordo, x)))
  
  # Remove groups without any descendants 
  ## In this case is 'Biological anomaly without phynotypic characterization'
  disease_with_no_descendants = sapply(disease_descendants, length)==0
  
  disease_descendants = disease_descendants[!disease_with_no_descendants]
  
  if(!return_genelist){
    return(disease_descendants)
  } else{
    terms_disease = terms_disease[!disease_with_no_descendants]
    
    # label the orphagroup with names
    orphaGroup_df = tibble(orphaGroupName = termLabel(terms_disease), orphaGroupID = termId(terms_disease))
    #disease_descendants_df = left_join(disease_descendants_df, orphaGroup_df)
    
    # a list containging disease associated genes
    disease_genes = lapply(disease_descendants, function(x){unique(unlist(combinedDiseaseGeneList[termId(x)]))
    })
    
    names(disease_genes) = termLabel(terms_disease)
    
    # wrap up as a data frame and wrote the tsv file with 
    # col1: disease name
    # col2: gene names, separated by semicolon ;
    #############################
    gene_disease_orpha = tibble(id = termId(terms_disease), name = termLabel(terms_disease), genes = disease_genes, n_genes = sapply(disease_genes, length))
    
    return(gene_disease_orpha)
  }
}

###############
# main: perform the retrieval of the associated genes given different roots
###############

#1: using 'Rare genetic disease' ( "Orphanet:98053")

gene_disease_orpha_genetic = retrieve_all_gene_association(root_term = "Orphanet:98053")
gene_disease_orpha_genetic = gene_disease_orpha_genetic %>% rowwise() %>% mutate(all_genes = paste(genes, collapse  = ";")) %>% select(name, all_genes)

write_tsv(gene_disease_orpha_genetic, "../data/table_disease_gene_assoc_orphanet_genetic.tsv")


#2: using 'group of disorders' (Orphanet:377794) - this is the uppermost term, consisting rare non-genetic and rare-genetic
gene_disease_orpha = retrieve_all_gene_association(root_term =  "Orphanet:377794")
gene_disease_orpha = gene_disease_orpha %>% rowwise() %>% mutate(all_genes = paste(genes, collapse  = ";")) %>% select(name, all_genes)

write_tsv(gene_disease_orpha, "../data/table_disease_gene_assoc_orphanet.tsv")





