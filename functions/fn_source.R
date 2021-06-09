##########################################
## Translational table from entrez to gene name
##########################################

# This function can take any of the columns(org.Hs.eg.db) as type and keys as long as the row names are in the format of the keys argument
IDconvert <- function(genelist, from="ENTREZID", to="SYMBOL"){
  #note from and to is taken from Affymatrix Human Genome U95 Set Annotation data; library('hgu95av2.db')
  # columns(hgu95av2.db) shows all possibility for to and from
  #  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  # [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  # [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  # [21] "PROBEID"      "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
  # [26] "UNIGENE"      "UNIPROT"   
  
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  
  geneSymbols <- mapIds(org.Hs.eg.db, keys=genelist, column=to, keytype=from, multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  #inds <- which(!is.na(geneSymbols))
  #found_genes <- geneSymbols[inds]
  
  return(geneSymbols)
}

# same as previous but lead with matrix row names
# This function can take any of the columns(org.Hs.eg.db) as type and keys as long as the row names are in the format of the keys argument
IDconvert_matrix <- function(df,  from="ENTREZID", to="SYMBOL"){
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  
  geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(df), column=to, keytype=from, multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}


############################
# dlapply: same as lapply but recursive so the function applies to all sublists (level 2)
############################
## define a function for applying on two levels of list
dlapply = function(list, func){
  nested_output = lapply(list, function(l1) 
    lapply(l1, function(l2) func(l2)))
  return(nested_output)
}



##############################
# create ggplot2 custom theme
###############################

theme_new <- function(base_size = 12, base_family = "Helvetica", xtilt = 0){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      #line = element_line(colour="black"),
      #text = element_text(colour="black"),
      axis.title = element_text(size = 14),
      #axis.text = element_text(colour="black", size=8),
      #strip.text = element_text(size=12),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.grid = element_blank(),   
      panel.border = element_rect(fill = NA, colour = "black", size=1),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}
