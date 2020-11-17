cluster_embeddings <- function(embs,file_info,method){
  library(tidyverse)
  library(clustree)
  library(Rtsne)
  library(ggplot2)
  library(umap)
}

prepape_embs <- function(emb, type = "unweighted", file_info, labels, keep_one, ave, n_emb){
  library(tidyverse)
  file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)
  # remove duplicate sig id embeddings
  id <- which(as.character(emb$emb) %in% as.character(file_info$emb))
  emb <- emb[id,]
  file_info <- file_info[which(as.character(file_info$emb) %in% as.character(emb$emb)),]
  # add sig id info to embs
  emb <- left_join(emb,file_info,by=c("emb"="emb"))
  # add label info to embs
  labels <- labels %>% group_by(moa_v1) %>% mutate(count = n_distinct(rdkit_graph)) %>% ungroup()
  labels <- labels %>% group_by(rdkit_graph) %>% filter(count == max(count)) %>% ungroup()
  labels <- labels %>% group_by(rdkit_graph) %>% sample_n(1) %>% ungroup()
  emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
  # keep one emb for each sig id
  if (keep_one) {
    emb <- emb %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
  }
  if (ave) {
    aver <- aggregate(emb[, 2:(n_emb+1)], list(emb$sig_id), mean)
    file_info_1 <- file_info %>% dplyr::select(sig_id,rdkit,cell_id) %>% unique()
    emb <- left_join(aver,file_info_1,by = c("Group.1"="sig_id"))
    emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
    colnames(emb)[1] <- "sig_id"
  }
  return(emb)
}

file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
labels <- readRDS("data/cmap/labels/labels_first_pass.rds")
id_topo <- which(grepl(pattern = "topoiso",x = labels$moa))
labels$moa_v1[id_topo] <- "topoisomerase_inhibitor"
# embs
library(data.table)
test <- fread("embeddings/graph2vec/emb_activity_1_epoch.csv")
colnames(test)[1] <- "emb"
n_emb = (ncol(test)-1)
emb_proc <- prepape_embs(emb = test,file_info = file_info,labels = labels,keep_one = F ,ave = T,n_emb = (ncol(test)-1))

# clustering kmeans

k <- seq(1,60,5)
clusters <- matrix(666,nrow = nrow(emb_proc),ncol = length(k))
colnames(clusters) <- as.character(k)
wss <- NULL
for (i in 1:length(k)) {
  test <- kmeans(scale(emb_proc[,2:(n_emb+1)]),centers = k[i],nstart = 50,iter.max = 40)
  clusters[,i] <- test$cluster
  colnames(clusters)[i] <- paste0("k",k[i])
  wss[i] <- test$tot.withinss
  print(i)
}

png(file="clustree.png",width=12,height=12,units = "in",res=300)
clustree(clusters, prefix = "k")
dev.off()


test <- kmeans(scale(emb_proc[,2:(n_emb+1)]),centers = 19,nstart = 50,iter.max = 40)
cluster <- test$cluster
emb_proc$cluster <- cluster
tsne_emb <- Rtsne(scale(emb_proc[,2:(n_emb+1)]), dims = 2, perplexity=85, verbose=TRUE, max_iter = 1000,initial_dims = 10,check_duplicates = F)

#fix moa, cell and disease area for vis

moa <- emb_proc %>%filter(!is.na(moa_v1)) %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count)) 
moa_vis <- moa$moa_v1[c(1,2,3,4,5,6,8,9,10)]
emb_proc$topmoa <- NA
emb_proc$topmoa[which(emb_proc$moa_v1 %in% moa_vis)] <- "moa"
emb_proc <- emb_proc %>% mutate(topmoa = if_else(condition = topmoa == "moa",true = moa_v1,false = topmoa))

cell <- emb_proc %>%filter(!is.na(cell_id)) %>% group_by(cell_id) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count)) 
cell_vis <- cell$cell_id[1:8]
emb_proc$cell_vis <- NA
emb_proc$cell_vis[which(emb_proc$cell_id %in% cell_vis)] <- "cell"
emb_proc <- emb_proc %>% mutate(cell_vis = if_else(condition = cell_vis == "cell",true = cell_id,false = as.factor(cell_vis)))

dis <- emb_proc %>%filter(!is.na(disease_area)) %>% group_by(disease_area) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count)) 
dis_vis <- dis$disease_area[2:9]
emb_proc$dis_vis <- NA
emb_proc$dis_vis <- as.factor(emb_proc$dis_vis)
emb_proc$dis_vis[which(emb_proc$disease_area %in% dis_vis)] <- "dis"
emb_proc <- emb_proc %>% mutate(dis_vis = if_else(condition = dis_vis == "dis",true = disease_area,false = as.character(dis_vis)))

df_emb <- data.frame(V1 = tsne_emb$Y[,1], V2 =tsne_emb$Y[,2], label = as.factor(emb_proc$topmoa),cluster = as.factor(emb_proc$cluster),
                     disease = as.factor(emb_proc$dis_vis),cell = as.factor(emb_proc$cell_vis))
df_emb$label <- factor(df_emb$label,levels = levels(df_emb$label),labels = c("ATP synthesis inhibitor",
                                                                             "CDK inhibitor","HDAC inhibitor","HSP inhibitor",
                                                                             "Mtor Inhibitor","NFKB inhibitor","PI3K inhibitor",
                                                                             "Protein synthesis inhibitor","Topoisomerase inhibitor"))
gtsne_clusters <- ggplot(df_emb, aes(V1, V2))+
  geom_point(aes(color = cluster),show.legend = T,size = 1) +
  geom_text(aes(label = cluster),size = 2)+
  xlab("X")+ylab("Y")+ 
  scale_color_discrete(name="Cluster")+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
gtsne_clusters



map <- umap(scale(emb_proc[,2:(n_emb+1)]))
df_map <- data.frame(V1 = map$layout[,1], V2 = map$layout[,2], label = as.factor(emb_proc$topmoa),cluster = as.factor(emb_proc$cluster))
df_map$label <- factor(df_map$label,levels = levels(df_map$label),labels = c("ATP synthesis inhibitor",
                                                                             "CDK inhibitor","HDAC inhibitor","HSP inhibitor",
                                                                             "Mtor Inhibitor","NFKB inhibitor","PI3K inhibitor",
                                                                             "Protein synthesis inhibitor","Topoisomerase inhibitor"))
umap_clusters <- ggplot(df_map, aes(V1, V2))+
  geom_point(aes(color = cluster),show.legend = T) +
  
  xlab("X")+ylab("Y")+ 
  scale_color_discrete(name="Cluster")+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
umap_clusters

umap_moa <- ggplot(df_map, aes(V1, V2))+
  geom_point(data= df_map %>% filter(is.na(label)),color = "lightgrey")+
  geom_point(data= df_map %>% filter(!is.na(label)),aes(color = label),show.legend = T) +
  xlab("X")+ylab("Y")+ 
  scale_color_manual(name="MoA",values = c("#984ea3","#f781bf","#ff7f00","#e41a1c","#4daf4a","#ffff33","#377eb8","#a65628","black"))+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
umap_moa+theme(legend.position = "bottom")

gtsne_moa <- ggplot(df_emb, aes(V1, V2))+
  geom_point(data= df_emb %>% filter(is.na(label)),color = "lightgrey",size = 1)+
  geom_point(data= df_emb %>% filter(!is.na(label)),aes(color = label),show.legend = T,size = 1) +
  xlab("X")+ylab("Y")+ 
  scale_color_manual(name="MoA",values = c("#984ea3","#f781bf","#ff7f00","#e41a1c","#4daf4a","#ffff33","#377eb8","#a65628","black"))+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
gtsne_moa+theme(legend.position = "bottom")
gtsne_moa+theme(legend.position = "none")

gtsne_cell <- ggplot(df_emb, aes(V1, V2))+
  geom_point(data= df_emb %>% filter(is.na(cell)),color = "lightgrey",size = 1)+
  geom_point(data= df_emb %>% filter(!is.na(cell)),aes(color = cell),show.legend = T,size = 1) +
  xlab("X")+ylab("Y")+ 
  scale_color_manual(name="Cell line",values = c("#984ea3","#f781bf","#ff7f00","#e41a1c","#4daf4a","#ffff33","#377eb8","#a65628"))+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
gtsne_cell+theme(legend.position = "bottom")

gtsne_dis <- ggplot(df_emb, aes(V1, V2))+
  geom_point(data= df_emb %>% filter(is.na(disease)),color = "lightgrey",size = 1)+
  geom_point(data= df_emb %>% filter(!is.na(disease)),aes(color = disease),show.legend = T,size = 1) +
  xlab("X")+ylab("Y")+ 
  scale_color_manual(name="Disease area",values = c("#984ea3","#f781bf","#ff7f00","#e41a1c","#4daf4a","#ffff33","#377eb8","#a65628"))+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)
gtsne_dis+theme(legend.position = "bottom")

png(file="disease_tsne.png",width=20,height=13,units = "cm",res=600)
gtsne_dis+theme(legend.position = "none")
dev.off()

png(file="clusters_tsne.png",width=20,height=20,units = "cm",res=600)
gtsne_clusters+theme(legend.position = "none")
dev.off()

png(file="moa_tsne.png",width=20,height=20,units = "cm",res=600)
gtsne_moa+theme(legend.position = "none")
dev.off()
