# Process the AUC from 10-fold cross-validation to include single most significant layer to a respective disease group
# Ize Buphamalai
# Jul 2021

#' Original AUC data include results on three sets of networks: PPI, all, signif
#' Inclusion of the most significant network (single layer benchmark)
#' While also update the 'signif' results - to incorporate the core co-expression layers


# Performance given by single most significant layer
AUC_folds_singleSignif <- readRDS("../cache/fold_cv_processed_results_revision_singleMostSignifLayer.RDS")

#AUC_folds_singleSignif <- AUC_folds_singleSignif %>% 
#  group_by(name, label) %>% 
#  summarise(mean = mean(AUC), med = median(AUC), Q25 = quantile(AUC, 0.25), Q75 = quantile(AUC, 0.75))

network_used <- AUC_folds_singleSignif %>% select(name, label) %>% distinct()

AUC_folds_singleSignif$label = "SingleMostSignif"

# Significant layers (update those with coex) - 

# compute median and 25% and 75% quantile
AUC_folds_signif_updated <- readRDS("../cache/fold_cv_processed_results_revision.RDS") %>% 
  filter(label %in% c("all", "signif_withCoexCore")) %>%
  mutate(label = factor(label,levels = c("all", "signif_withCoexCore"), 
                        labels = c("All layers", "Significant layers")))

disease_withUpdatedData <- AUC_folds_signif_updated %>% filter(label == "Significant layers") %>% pull(name) %>% unique

# the signif network inly computed for those that coex_core does not change   
AUC_folds_signif_orig <- readRDS("../cache/fold_cv_processed_results.RDS") %>% filter((label=="significant networks" & !name %in% disease_withUpdatedData) | label == "PPI only") %>%
  mutate(label = factor(label, levels = c("significant networks", "PPI only"), labels = c("Significant layers", "PPI"))) %>% select(-name_short)


AUC_folds_summary_updated <- rbind(AUC_folds_signif_updated, AUC_folds_singleSignif) %>% 
  group_by(name, label) %>% 
  summarise(mean = mean(AUC), med = median(AUC), Q25 = quantile(AUC, 0.25), Q75 = quantile(AUC, 0.75)) #%>%
#  mutate(name_short = factor(name, levels = levels(name), labels = short_names))

AUC_folds_summary_updated <- rbind(AUC_folds_summary_updated,AUC_folds_signif_orig)

# reorder based on performance
order_performance <- AUC_folds_summary_updated %>% filter(label == "Significant layers") %>% arrange(mean) %>% pull(name)


AUC_folds_summary_updated_plot <- AUC_folds_summary_updated %>%
  mutate(name = factor(name, levels = order_performance),
         name = fct_recode(name, `Rare genetic gynecological diseases`="Rare genetic gynecological and obstetrical diseases", 
                           `Rare genetic developmental defect..`="Rare genetic developmental defect during embryogenesis",
                           `Rare genetic rheumatologic disease`="Rare genetic systemic or rheumatologic disease"))

saveRDS(AUC_folds_summary_updated_plot, "../cache/fold_cv_processed_results_revisioned_processed.RDS")
