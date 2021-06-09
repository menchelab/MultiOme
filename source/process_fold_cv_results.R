# Process results from LOOCV cross validation
# Ize Buphamalai


pacman::p_load("cvAUC")

## #################
## load the data

# result from using significant networks
rank_df_all_folds <- readRDS("../cache/LCC_CV_ranking_rare_genetic_diseases_no_LCC_recompute.RDS")
rank_df_all_folds <- lapply(rank_df_all_folds, function(x) lapply(x, function(df) df[, c("trueset", "RankArithmP")]))

# result from using PPI alone
rank_df_all_folds_only_PPI <- readRDS("../cache/LCC_CV_ranking_rare_genetic_diseases_only_PPI.RDS")
# Since using one network there is one ranking column
rank_df_all_folds_only_PPI <- lapply(rank_df_all_folds_only_PPI, function(x) lapply(x, function(df) df[, c("trueset", "rank")]))

# results from all networks
rank_df_all_folds_allnets <-readRDS("../cache/LCC_CV_ranking_rare_genetic_diseases_allnets.RDS")
rank_df_all_folds_allnets_used <- lapply(rank_df_all_folds_allnets, function(x) lapply(x, function(df) df[, c("trueset", "RankArithmP")]))


############
# Process the AUC values
# function to compute AUC from rank

auc_output<-function(rank_df){
  #' @input rank_df: a two column data frame, with label and prediction respectively
  rank_label<-lapply(rank_df, function(x) x[,1])
  rank_prediction<-lapply(rank_df, function(x) -x[,2])
  
  out <- cvAUC(rank_prediction, rank_label)
  return(out)
}


## compute AUC from the rank through CVAUC - combining cross validation results

aucvals_signif <-lapply(rank_df_all_folds, auc_output)
aucvals_PPI <-lapply(rank_df_all_folds_only_PPI, auc_output)
aucvals_all <-lapply(rank_df_all_folds_allnets_used, auc_output)

all_diseases <- names(aucvals_all)


##################################
# AUC plot per disease, compared each method
# uncomment this section to make plots

# all diseases in one plot

# palettes <- c('#66c2a5','#fc8d62','#8da0cb')
## signif only
# pdf(paste0("../Figs/cvROC/all_rare_genetic_diseases_no_LCC_remeasure.pdf"), width=12, height=10)
# par(mfrow=c(5,6))
# {
#   for(i in all_diseases)  #Plot fold AUCs
#   { plot(aucvals_signif[[i]]$perf, col="grey82", lty=3, main=paste(i))
#     
#     #Plot CV AUC
#     plot(aucvals_signif[[i]]$perf, col="red", avg="vertical", add=TRUE)}
# }
# dev.off()
# 
# 
# ## PPI only
# pdf(paste0("../Figs/cvROC/all_rare_genetic_diseases_PPI_only.pdf"), width=12, height=10)
# par(mfrow=c(5,6))
# {
#   for(i in all_diseases)  #Plot fold AUCs
#   { plot(aucvals_PPI[[i]]$perf, col="grey82", lty=3, main=paste(i))
#     
#     #Plot CV AUC
#     plot(aucvals_PPI[[i]]$perf, col="red", avg="vertical", add=TRUE)}
# }
# dev.off()
# 
# ## allnets
# pdf(paste0("../Figs/cvROC/all_rare_genetic_diseases_all_networks.pdf"), width=12, height=10)
# par(mfrow=c(5,6))
# {
#   for(i in all_diseases)  #Plot fold AUCs
#   { plot(aucvals_all[[i]]$perf, col="grey82", lty=3, main=paste(i))
#     
#     #Plot CV AUC
#     plot(aucvals_all[[i]]$perf, col="red", avg="vertical", add=TRUE)}
# }
# dev.off()

## AUC plot

# pdf(paste0("../Figs/cvROC/all_rare_genetic_diseases_compared_all_methods.pdf"), paper = "a4")
# par(mfrow=c(7,4))
# {
#   for(i in all_diseases)  #Plot fold AUCs
#   { 
#     # signif
#     plot(aucvals_signif[[i]]$perf, col=palettes[1],  avg="vertical", main=str_remove_all(i, "Rare |genetic | disease"))
#     
#     # all networks
#     plot(aucvals_all[[i]]$perf, col=palettes[2], avg="vertical", add=TRUE)
#     
#     # PPI alone
#     plot(aucvals_PPI[[i]]$perf, col=palettes[3], avg="vertical", add=TRUE)
#     
#   }
# }
# dev.off()
# 
# 
# 
# for(i in all_diseases)  #Plot fold AUCs
# { 
#   pdf(paste0("../Figs/cvROC/by_disease/",i,".pdf"), width = 4, height = 4)
#   par(mar = c(2,2,2,1)+0.1)
#   # signif
#   plot(aucvals_signif[[i]]$perf, col=palettes[3],  avg="vertical", main=str_remove_all(i, "Rare |genetic | disease| diseases"), xlab = "", ylab ="", font.main = 8, lwd = 3)
#   
#   # all networks
#   plot(aucvals_all[[i]]$perf, col=palettes[1], avg="vertical", add=TRUE, lwd = 2)
#   
#   # PPI alone
#   plot(aucvals_PPI[[i]]$perf, col=palettes[2], avg="vertical", add=TRUE, lwd = 2)
#   
#   dev.off()
# }


######
# process auc data into data frame for post analysis
# AUC plot for all diseases
auc_to_df <- function(auc_object, label = NULL){
  AUC_folds<-lapply(auc_object, function(x) x$fold.AUC)
  
  AUC_folds<- reshape2::melt(AUC_folds)
  colnames(AUC_folds)<-c("AUC", "name")
  
  AUC_folds<-AUC_folds %>% mutate(name=fct_reorder(name, AUC))
  
  if(!is.null(label)){
    AUC_folds$label <- label
  }
  return(AUC_folds)
}


### AUC plot for all diseases and all methods
AUC_folds_signif <- auc_to_df(aucvals_signif, label = "significant networks")
AUC_folds_allnets <- auc_to_df(aucvals_all, label = "all networks")
AUC_folds_PPI <- auc_to_df(aucvals_PPI, label = "PPI only")

AUC_folds <- bind_rows(AUC_folds_signif, AUC_folds_allnets, AUC_folds_PPI)

# add short names for easy labelling
short_names = c(
  "Rare genetic developmental defect...", 
  "Rare genetic endocrine disease", 
  "Rare genetic neurological disorder", 
  "Rare chromosomal anomaly", 
  "Rare inborn errors of metabolism", 
  "Rare genetic eye disease", 
  "Rare genetic bone disease", 
  "Rare genetic tumor", 
  "Rare genetic immune disease", 
  "Rare genetic skin disease", 
  "Genetic otorhinolaryngologic disease", 
  "Rare genetic cardiac disease", 
  "Rare genetic renal disease", 
  "Genetic infertility", 
  "Rare genetic urogenital disease", 
  "Rare genetic gastroenterological disease", 
  "Rare genetic hepatic disease", 
  "Rare genetic gynecological diseases", 
  "Rare genetic hematologic disease", 
  "Rare genetic rheumatologic disease", 
  "Ciliopathy", 
  "Inherited cancer-predisposing syndrome", 
  "Rare genetic odontologic disease", 
  "Laminopathy", 
  "Rare genetic vascular disease", 
  "RASopathy")



# compute median and 25% and 75% quantile
AUC_folds_summary <- AUC_folds %>% 
  group_by(name, label) %>% 
  summarise(mean = mean(AUC), med = median(AUC), Q25 = quantile(AUC, 0.25), Q75 = quantile(AUC, 0.75)) %>%
  mutate(name_short = factor(name, levels = levels(name), labels = short_names))
  


saveRDS(AUC_folds_summary, "../cache/fold_cv_processed_results.RDS")
