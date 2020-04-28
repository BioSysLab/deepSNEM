library(tidyverse)

### cmap drugs

cmap <- readRDS("data/cmap/util_files/pert_id_to_rdkit.rds")


# read broad repurposing data

broad_repo <- readRDS("data/cmap/repo_hub/broad_repurposing.rds")



#### read the quality 1 smiles for which we hve the graphs

data_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
data <- readRDS("data/graph_info_df/file_info_nodups.rds")

data <- data %>% select(rdkit) %>% filter(!is.na(rdkit)) %>% unique()
data_dups <- left_join(data_dups,cmap, by = "pert_id")
data_dups <- data_dups %>% select(rdkit) %>% filter(!is.na(rdkit)) %>% unique()

### rdkit with graph available
q1_rdkits <- bind_rows(data,data_dups) %>% unique()

sims <- read.csv("data/cmap/repo_hub/deepsnem_graph_repo_sims.csv") 
sims <- sims[,-1]
sims <- sims > 0.99
sims <- sims+0
rownames(sims) <- as.character(q1_rdkits$rdkit)
sims <- sims[which(rowSums(sims)!=0),]

# add labels to rdkit wiith graphs
labels <- data.frame(matrix(0,nrow=nrow(sims),ncol=ncol(broad_repo)+1))
colnames(labels) <- c("rdkit_graph","rdkit_broad","moa","target","disease_area","indication")
broad_repo$rdkit <- as.character(broad_repo$rdkit)
broad_repo$moa <- as.character(broad_repo$moa)
broad_repo$target <- as.character(broad_repo$target)
broad_repo$disease_area <- as.character(broad_repo$disease_area)
broad_repo$indication <- as.character(broad_repo$indication)
for (i in 1:nrow(sims)) {
  id <- which(sims[i,] == 1)[1]
  labels[i,"rdkit_graph"] <- rownames(sims)[i]
  labels[i,2:ncol(labels)] <- as.character(broad_repo[id,])
}

saveRDS(labels,"data/cmap/labels.rds")
write.csv(labels,"data/cmap/labels.csv")
moa <- labels %>% group_by(moa) %>% summarise(count = n()) %>% arrange(desc(count))
target <- labels %>% group_by(target) %>% summarise(count = n()) %>% arrange(desc(count))
disease_area <- labels %>% group_by(disease_area) %>% summarise(count = n()) %>% arrange(desc(count))
indication <- labels %>% group_by(indication) %>% summarise(count = n()) %>% arrange(desc(count))
