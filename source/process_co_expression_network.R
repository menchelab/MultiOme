################################################
#' the disparity filter network with specificity
#' set a cutoff of edges after applying disparity filter
# rationale: edges strongly identified to be in many networks are perhaps interactions between house keeping genes and not a part of tissue specific mechanisms



# 0 - load the co-expression data --------------
coex_folder = "./networks/gtex_coexpresssion/edgelists/coex_dispfilter_001_cor075/"

dat = list()
for(i in list.files(coex_folder)){
  tissue_name = str_remove(i, ".tsv")
  
  dat[[tissue_name]] = read_tsv(paste0(coex_folder, i), col_types = "cc")

}


# 1 - combine edges from all tissues and count occurrences --------------
dat_combn = bind_rows(dat)
dat_combn= dat_combn %>% mutate_if(is.character,as.factor)
dat_combn_unique = dat_combn %>% group_by(A,B) %>% count()

saveRDS(dat_combn_unique, "./cache/coexpression_raw_edge_counts.RDS")


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


##############

# 
# # new coex by group (without disparity filter)
# 
# coex_folder = "~/Documents/projects/Multiome/networks/gtex_coexpresssion/edgelists/new_coex_by_group//"
# 
# dat = list()
# for(i in list.files(coex_folder)){
#   tissue_name = str_remove(i, ".tsv")
#   
#   dat[[tissue_name]] = read_tsv(paste0(coex_folder, i), col_types = "cc")
#   colnames(dat[[tissue_name]]) = c("A","B")
# }
# 
# dat_combn = bind_rows(dat)
# dat_combn= dat_combn %>% mutate_if(is.character,as.factor)
# colnames(dat_combn) =  c("A","B")
# dat_combn_unique = dat_combn %>% group_by(A,B) %>% count()
# 
# # set a cut off for a maixum network to tolerate
# cutoff = 5
# dat_combn_pass = dat_combn_unique[dat_combn_unique$n <= cutoff, c("A","B")]
# 
# # publish the new sets of networks
# for(i in names(dat)){
#   print(i)
#   dat_pass = inner_join(dat[[i]], dat_combn_pass,  by = c("A", "B"))
#   write_tsv(dat_pass, paste0("~/Documents/projects/Multiome/networks/gtex_coexpresssion/edgelists/coex_group38_cor075_cut5/",i,".tsv"))
# }
# 
# ######## cutoff of ten
# 
# # set a cut off for a maixum network to tolerate
# cutoff = 10
# dat_combn_pass = dat_combn_unique[dat_combn_unique$n <= cutoff, c("A","B")]
# 
# # publish the new sets of networks
# for(i in names(dat)[34:38]){
#   print(i)
#   dat_pass = inner_join(dat[[i]], dat_combn_pass,  by = c("A", "B"))
#   write_tsv(dat_pass, paste0("~/Documents/projects/Multiome/networks/gtex_coexpresssion/edgelists/coex_group38_cor075_cut5/",i,".tsv"))
# }
