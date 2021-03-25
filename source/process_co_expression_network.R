################################################
#' coe-expression network post-processing characterisation
#' 
#' The co-expression edges, after disparity filter were labelled with 
#' 1) whether they are from essential gene pairs
#' 2) how many tissues they appear in
#' 3) set a cutoff of for specificity: edges appear for over a certain number of tissues does not show tissue specific signal
#' this number was set to 5 (allowing some redundancy for tissues that are related)


# 0 - load the co-expression data --------------
coex_folder = "./networks/gtex_coexpresssion/edgelists/coex_dispfilter_001_cor075/"

dat = list()
for(i in list.files(coex_folder)){
  tissue_name = str_remove(i, ".tsv")
  
  dat[[tissue_name]] = read_tsv(paste0(coex_folder, i), col_types = "cc")

}


# 1 - combine edges from all tissues and count occurrences --------------
coex_el_sum = bind_rows(dat)
coex_el_sum = coex_el_sum %>% mutate_if(is.character,as.factor)
# count the number of tissue each edge emerges
coex_el_sum = coex_el_sum %>% group_by(A,B) %>% count()


#2 - incorporate essential genes data
source("../functions/fn_source.R")
essential_genes = read_tsv("../data/OGEE_esential_genes_20190416.txt")
essential_genes$geneSymbol = IDconvert(essential_genes$locus, from = "ENSEMBL", to = "SYMBOL")
essential_genes_sym = essential_genes %>% filter(essential == "E") %>% pull(geneSymbol) %>% unique()

# annodate whether each node is an essential gene, and whether the pair is both essential
coex_el_sum$A.essential = coex_el_sum$A %in% essential_genes_sym
coex_el_sum$B.essential = coex_el_sum$B %in% essential_genes_sym
coex_el_sum$essential_edge_score = coex_el_sum$A.essential + coex_el_sum$B.essential

# group occurrences by increment of five tissues
coex_el_sum$n_binned = cut(coex_el_sum$n, breaks = c(seq(0,30,5), 38))

saveRDS(coex_el_sum, "../cache/coexpression_raw_edge_counts.RDS")

#- store the value by group

# change label from domain-range expression for readibility
changelabel = function(x){gsub("\\(","", gsub("\\]", "", gsub(",","-", x)))} 

coex_el_sum_grouped = coex_el_sum  %>% ungroup %>% group_by(n_binned, essential_edge_score) %>% count()

coex_el_sum_grouped = coex_el_sum_grouped %>% mutate(
  n_binned_relabel = factor(changelabel(n_binned), levels =  changelabel(levels(n_binned))),
  score = as.character(essential_edge_score))

saveRDS(coex_el_sum_grouped, "../cache/coexpression_edge_counts_by_group.RDS")


# 2 - set a cut off for a maixum network to tolerate --------------

### 2.1 - allowing max 5 networks 
cutoff = 5
dat_combn_pass = dat_combn_unique[dat_combn_unique$n <= cutoff, c("A","B")]

# publish the new sets of networks
for(i in names(dat)){
  print(i)
  dat_pass = inner_join(dat[[i]], dat_combn_pass,  by = c("A", "B"))
  write_tsv(dat_pass, paste0("./networks/gtex_coexpresssion/edgelists/coex_dispfilter_001_cor075_cut5/",i,".tsv"))
}


### 2.2 - allowing max 10 networks 
cutoff = 10
dat_combn_pass = dat_combn_unique[dat_combn_unique$n <= cutoff, c("A","B")]

# publish the new sets of networks
for(i in names(dat)){
  print(i)
  dat_pass = inner_join(dat[[i]], dat_combn_pass,  by = c("A", "B"))
  write_tsv(dat_pass, paste0("./networks/gtex_coexpresssion/edgelists/coex_dispfilter_001_cor075_cut10/",i,".tsv"))
}


