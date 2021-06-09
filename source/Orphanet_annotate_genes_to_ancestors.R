# Associate disease genes to ancestor terms
# Ize Buphamalai
# Updated: 09.04.2020


# input: ontology structure: ordo (EBI lookup service rols object)
# input: gene-orphanet disease association

# The code find all ancestors for all terms and sequentially associate genes to upper levels on ontology

#################
## Load required packages and data
################
library(pacman)
pacman::p_load("rols", "stringr")

#Read orphanet gene association data (immediate association)
source("./read_orphanet_gene_association_data.R")


#################
## Ontological data extraction
################
# extract the ontology structure from Ontology Lookup service
print("read the ontology from the EBI ontology lookup service")

ordo = Ontology("ordo")
terms = terms(ordo)

# store all term IDs
all_terms = termId(terms)

terms_df = as(terms, "data.frame")

## translational table between diseases with genes and disease groups
terms_with_genes = terms[termId(terms) %in% orphaGeneDisease$orphaID]

##############
# Modify a function to collect the ancestor terms from the ontology
# is_a: natural recognition of ancestor search
# part_of: a disease is a part of a disease group, but not directly detected from ancestor search on rols package. 
# the following part of the code is to search two levels up on 'part of' branch
###############

# create a new function based on PartOf
# suppress messsage for terms with no associated 'PartOf' terms

# following is the copy of PartOf from the rols package, modified by removing the message line
# and replace it to the function

partOf_suppress = function (id) 
{
  stopifnot(inherits(id, "Term"))
  url <- id@links$part_of[[1]]
  if (is.null(url)) {
    return(NULL)
  }
  x <- httr::GET(url)
  stop_for_status(x)
  cx <- content(x)
  ans <- lapply(cx[["_embedded"]][[1]], makeTerm)
  names(ans) <- sapply(ans, termId)
  Terms(x = ans)
}

assignInNamespace("partOf", partOf_suppress, ns="rols")
#fixInNamespace("partOf", "rols")

partOf_extend = function(term){
  firstLevel = partOf(term)
  secondLevel = lapply(firstLevel, function(x) partOf(x))
  results = c(unlist(firstLevel), unlist(secondLevel))
  results = unique(unlist(lapply(results, termId)))
  return(results)
}

# annotate 'partof' terms to the terms with associated genes
############
print(" annotate 'partof' terms to the terms with associated genes")

# compute ' part of' terms - results are terms that the disease is a 'part of'
terms_with_genes_parts_extend = lapply(terms_with_genes,  function(x) suppressMessages(partOf_extend(x)))

# remove terms labelled as NULL
terms_with_genes_parts_extend_rm_null = terms_with_genes_parts_extend[!sapply(terms_with_genes_parts_extend, is.null)]


# create a function to retrieve all ancestors of a term
###############
print("retrieve all ancestors of all terms")

all_ancestors = function(termlist){
  if(is.character(termlist)){
    termlist = terms[all_terms %in% termlist]
  }
  IDancestors = lapply(termlist, function(y) termId(ancestors((y))))
  IDancestors = unique(unlist(IDancestors))
  return(IDancestors)
}

# apply the function -- this process might take very long 
terms_with_genes_ancestors = lapply(terms_with_genes_parts_extend_rm_null,  all_ancestors)


#############
## get the parents info

all_children = function(termlist){
  if(is.character(termlist)){
    termlist = terms[all_terms %in% termlist]
  }
  IDchildren = lapply(termlist, function(y) termId(children((y))))
  IDchildren = unique(unlist(IDchildren))
  return(IDchildren)
}
# apply the function -- this process might take very long 

# descendants of the rare genetic diseases
termlist = descendants(terms[all_terms == "Orphanet:98053"][[1]])

terms_children = lapply(termlist, function(x) suppressMessages(children(x)))


term_childern_chr = terms_children[!sapply(terms_children, length) == 0]
  
term_childern_chr =  lapply(term_childern_chr, termId)
term_childern_df = melt(term_childern_chr)
colnames(term_childern_df) = c("to", "from")


disease_descendants = lapply(terms[all_terms == "Orphanet:98053"], descendants)

#################
# store the list of genes associated to disease group
#################
print("annotate disease associate genes cumulatively to ancestors")

# Expanded disease Gene List: associate gene to ancestor terms as well
expandedDiseaseGeneList = list()

for(i in names(orphaGeneDiseaseList)){
  extended_terms = terms_with_genes_parts_extend_rm_null[[i]]
  if(!is.null(extended_terms)){
    #extended_terms = termId(extended_terms)
    for(j in extended_terms){
      if(!is_empty(expandedDiseaseGeneList[[j]])){
        expandedDiseaseGeneList[[j]] = append(expandedDiseaseGeneList[[j]], orphaGeneDiseaseList[[i]])
      }
      else{
        expandedDiseaseGeneList[[j]] = orphaGeneDiseaseList[[i]]
      }
    }
  }
}

# combine expanded association with the original association
combinedDiseaseGeneList = c(expandedDiseaseGeneList, orphaGeneDiseaseList)

print("saving the results")

saveRDS(combinedDiseaseGeneList, file = "../cache/orphanet_disease_gene_association_with_ancestors.RDS")

##################


#################
# orphaID to name
#################
orphaID_to_name_df = read_csv("../data/orphanumbers_to_name.csv", col_names = c("orphaID", "orphaName"))
orphaID_to_name_df$orphaID = paste0("Orphanet:", orphaID_to_name_df$orphaID)

orpha_to_gene_df = melt(orphaGeneDiseaseList)
colnames(orpha_to_gene_df) = c("gene", "orphaID")

############################
# extract the children of 'rare genetic diseases' and the labels

rare_genetic_subterms <- termId(children(terms[all_terms == "Orphanet:98053"][[1]]))
rare_genetic_subterms_label <- sapply(rare_genetic_subterms, function(x) terms[[x]]@label)

# take only ancestors that are in rare genetic diseases
terms_with_genes_ancestors_in_rare_genetic <- lapply(terms_with_genes_ancestors, function(x) intersect(x, rare_genetic_subterms))
terms_with_genes_ancestors_in_rare_genetic_label <- lapply(terms_with_genes_ancestors_in_rare_genetic, function(x) paste(rare_genetic_subterms_label[x],collapse = ";"))

terms_rare_genetic_label_df <- melt(terms_with_genes_ancestors_in_rare_genetic_label)
colnames(terms_rare_genetic_label_df) <- c("roots", "orphaID")

#load orphaGeneDisease
source("./read_orphanet_gene_association_data.R")

orphaGeneDisease_with_ancestor <- left_join(orphaGeneDisease, terms_rare_genetic_label_df)
orphaGeneDisease_with_ancestor$genes <- sapply(orphaGeneDisease_with_ancestor$genes, function(x) paste(x, collapse = ";"))
orphaGeneDisease_with_ancestor <- left_join(orphaGeneDisease_with_ancestor, terms_df[,c("id", "label","description")], by = c("orphaID" = "id"))

write_tsv(orphaGeneDisease_with_ancestor[,c("orphaID","label","genes","n_genes", "description", "roots")], path = "../data/orphaNet_disease_gene_association_with_roots.tsv")

###################
## add the onset
#################
## Onset information
#################
onsetfiles = list.files("../data/orphanet_onset/")

# get orphanet IDs that the onset was mentioned
onset_id = list()

for(i in onsetfiles){
  onset_disease = readLines(paste0("../data/orphanet_onset/", i))
  onset_id[[i]] = terms_df %>% dplyr::filter(label %in% onset_disease) %>% pull(id)
}

orpha_to_onset_df = melt(onset_id); colnames(orpha_to_onset_df) = c("orphaID", "onset")


#################
## merge orpha - gene - onset

# join the data
orpha_gene_onset_df = full_join(orpha_to_onset_df, orpha_to_gene_df)

# make onset a ordered factor
orpha_gene_onset_df$onset = factor(orpha_gene_onset_df$onset, levels = c("antenatal", "neonatal",  "infancy", "childhood",  
                                                                         "adolescent", "adult", "all_ages", "elderly"), 
                                   ordered = T)

saveRDS(orpha_gene_onset_df, "../cache/orpha_gene_onset_df.RDS")


# 1.summarise for each disease
orpha_by_earliest_onset_df = orpha_gene_onset_df %>% 
  dplyr::filter(!is.na(onset), !is.na(gene)) %>%
  group_by(orphaID) %>% summarise(earliest_onset = min(onset))

# map the orphaName to it
orpha_by_earliest_onset_df = left_join(orpha_by_earliest_onset_df, orphaID_to_name_df) 


# 2. summarise based on earliest onset for each gene
gene_by_earliest_onset_df = orpha_gene_onset_df %>% 
  dplyr::filter(!is.na(onset), !is.na(gene)) %>%
  group_by(gene) %>% summarise(earliest_onset = min(onset))

# number of genes associated to each onset type
gene_by_earliest_onset_df %>% group_by(earliest_onset) %>% count()

# save the information
saveRDS(gene_by_earliest_onset_df, "../cache/gene_by_earliest_onset_orphanet.RDS")



