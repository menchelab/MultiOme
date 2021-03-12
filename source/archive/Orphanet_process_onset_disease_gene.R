# Gene disease and onset data
# Pisanu Ize Buphamalai
# Updated: 09.04.2019


#################
## Gene-disease associaton
#################
source("./read_orphanet_gene_association_data.R")


#################
# orphaID to name
#################
orphaID_to_name_df = read_csv("../data/orphanumbers_to_name.csv", col_names = c("orphaID", "orphaName"))
orphaID_to_name_df$orphaID = paste0("Orphanet:", orphaID_to_name_df$orphaID)

orpha_to_gene_df = melt(orphaGeneDiseaseList)
colnames(orpha_to_gene_df) = c("gene", "orphaID")

#################
## Onset information
#################
onsetfiles = list.files("../data/orphanet_onset/")

# get orphanet IDs that the onset was mentioned
onset_id = list()

for(i in onsetfiles){
  onset_disease = readLines(paste0("../data/orphanet_onset/", i))
  onset_id[[i]] = terms_df %>% filter(label %in% onset_disease) %>% pull(id)
}

orpha_to_onset_df = melt(onset_id); colnames(orpha_to_onset_df) = c("orphaID", "onset")


#################
## merge orpha - gene - onset
#################

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

