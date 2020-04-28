cluster_embeddings <- function(embs,file_info,method){
  library(tidyverse)
  library(clustree)
}

# file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)

# embs

emb <- read.csv("embeddings/graph2vec/emb_clustered_norm_500.csv")
colnames(emb)[1] <- "emb"
#emb <- emb %>% mutate(emb = str_remove_all(string = emb,pattern = ".csv"))

# keep only 1 emb per sig id for clustering
# keep the respective embeddings

emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
file_info <- file_info[which(as.character(file_info$emb) %in% as.character(emb$emb)),]

file_info <- file_info %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]

k <- seq(1,60)
clusters <- matrix(666,nrow = nrow(emb),ncol = length(k))
colnames(clusters) <- as.character(k)
wss <- NULL
for (i in 1:length(k)) {
  test <- kmeans(scale(emb[,-1]),centers = k[i],nstart = 50,iter.max = 40)
  clusters[,i] <- test$cluster
  colnames(clusters)[i] <- paste0("k",k[i])
  wss[i] <- test$tot.withinss
}
png(file="clustree.png",width=12,height=12,units = "in",res=300)
clustree(clusters, prefix = "k")
dev.off()
