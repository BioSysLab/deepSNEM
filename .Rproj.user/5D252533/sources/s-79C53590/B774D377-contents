prepape_embs <- function(emb, type = "unweighted", file_info, labels, keep_one, ave){
  library(tidyverse)
  file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)
  # remove duplicate sig id embeddings
  emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
  file_info <- file_info[which(as.character(file_info$emb) %in% as.character(emb$emb)),]
  # add sig id info to embs
  emb <- left_join(emb,file_info,by=c("emb"="emb"))
  # add label info to embs
  emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
  # keep one emb for each sig id
  if (keep_one) {
    emb <- emb %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
  }
  if (ave) {
    aver <- aggregate(emb[, 2:513], list(emb$sig_id), mean)
    file_info_1 <- file_info %>% select(sig_id,rdkit) %>% unique()
    emb <- left_join(aver,file_info_1,by = c("Group.1"="sig_id"))
    emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
    colnames(emb)[1] <- "sig_id"
  }
  return(emb)
}
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
emb <- read.csv("embeddings/ged_distance/ged_embs512_seen_4ep.csv")
emb <- emb[,-1]
colnames(emb)[1] <- "emb"
saveRDS(emb,"embeddings/ged_distance/ged_embs512_seen_4ep.rds")
labels <- readRDS("data/cmap/labels_first_pass.rds")
emb_proc <- prepape_embs(emb = emb,file_info = file_info,labels = labels,keep_one = F ,ave = T)


visualize_moa <- function(emb,output_dir,moa_n,emb_size,perpl_emb,iter,init_dim,name){
  library(tidyverse)
  library(Rtsne)
  addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_text(size = textSize), 
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }
  moa <- emb %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
  moa_vis <- moa$moa_v1[1:moa_n]
  emb <- emb[which(as.character(emb$moa_v1) %in% moa_vis),]
  tsne_emb <- Rtsne(scale(emb[,2:(emb_size+1)]), dims = 2, perplexity=perpl_emb, verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
  df_emb <- data.frame(V1 = tsne_emb$Y[,1], V2 =tsne_emb$Y[,2], label = as.factor(emb$moa_v1))
  gtsne <- ggplot(df_emb, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_discrete()
  png(file=paste0(output_dir,"/",as.character(name),"_moa_tsne.png"),width=9,height=9,units = "in",res=300)
  print(addSmallLegend(gtsne))
  dev.off()
}

visualize_moa(emb_proc,output_dir = "vis",moa_n = 10,emb_size = 512,
              perpl_emb = 10,init_dim = 80,iter = 1000,name = "test_ave_10perpl_1000_iter_10moa")
