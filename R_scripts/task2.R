emb_task2 <- function(test, method_name, file_info, distance_type, file_info_dups){
  # test is a dataframe where there is a column emb which is in the
  # form of the column emb of file info
  # file info, dataframe 
  # distance, character string, either "cosine" or "euclidian"
  
  # Task 2 validation embeddings
  library(tidyverse)
  library(lsa)
  library(Rtsne)
  file_info <- file_info %>% dplyr::select(-c(distil_id,distil_cc_q75,distil_nsample,distil_ss,median_drug_ranks,ngenes_modulated_up_lm,
                                              ngenes_modulated_dn_lm,tas,pct_self_rank_q25,is_exemplar,myc_ranks,counts,count.y))
  ids_keep <- which(file_info$emb %in% test$emb)
  file_info <- file_info[ids_keep,]
  sig_counts <- file_info %>% group_by(sig_id) %>% summarize(count_per_sig = n())
  file_info <- left_join(file_info,sig_counts) 
  file_info_dups$identifier <- as.character(file_info_dups$identifier)
  
  if (distance_type == "cosine") {
    # cosine is normalized to 0-1 universally
    distance_function <- function(df){
      cos <- cosine(t(df[,-1]),y = NULL)
      cos <- (cos + 1)/(2)
      cos_dist <- 1-cos
      cos_dist[lower.tri(cos_dist,diag = T)] <- 666
      rownames(cos_dist) <- as.character(df$emb)
      colnames(cos_dist) <- as.character(df$emb)
      cos_dist <- reshape2::melt(cos_dist)
      cos_dist <- cos_dist %>% filter(value != 666)
      return(cos_dist)
    }
  }
  
  if (distance_type == "euclidian") {
    #this is unnormalized
    distance_function <- function(df){
      eu_dist <- dist(df[,-1])
      eu_dist <- as.matrix(eu_dist)
      eu_dist[lower.tri(eu_dist,diag = T)] <- 666
      rownames(eu_dist) <- as.character(df$emb)
      colnames(eu_dist) <- as.character(df$emb)
      eu_dist <- reshape2::melt(eu_dist)
      eu_dist <- eu_dist %>% filter(value != 666)
      return(eu_dist)
    }
  }
  
  # Task 2 validation embeddings
  
  identifier <- unique(as.character(file_info_dups$identifier))
  dup_distances <- NULL
  for (i in 1:length(identifier)) {
    filt <- file_info_dups[which(file_info_dups$identifier %in% identifier[i]),]
    sigs_iden <- unique(as.character(filt$sig_id))
    
    for (j in 1:10) {
      dup_embs <- NULL
      for (k in 1:length(sigs_iden)) {
        filt2 <- filt %>% filter(sig_id == sigs_iden[k])
        filt_embs <- test[which(test$emb %in% filt2$emb),]
        dup_emb <- sample_n(filt_embs,1)
        dup_embs <- rbind(dup_emb,dup_embs)
      }
      if (nrow(dup_embs)>1) {
        dup_distance <- distance_function(dup_embs)
        dup_distances <- rbind(dup_distance,dup_distances)
      }
      
    }
  }
  
  sigs <- unique(as.character(file_info$sig_id))
  set.seed(1337)
  sigs_random <- sample(sigs,90)
  
  random_embs <- NULL
  for (i in 1:length(sigs_random)) {
    filt <- file_info %>% filter(sig_id == sigs_random[i])
    random_emb <- sample_n(test[which(test$emb %in% filt$emb),],1)
    random_embs <- rbind(random_emb,random_embs)
  }
  random_distances <- distance_function(random_embs)
  
  dup_distances$method <- method_name
  dup_distances$type <- "Duplicate signatures"
  
  random_distances$method <- method_name
  random_distances$type <- "Random signatures"
  
  task2 <- bind_rows(dup_distances,random_distances)
  task2$value <- task2$value/max(task2$value)
  return(task2)
  
}
library(data.table)
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
# dup file
file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
test <- fread("embeddings/graph2vec/emb_activity_1_epoch.csv")
colnames(test)[1] <- "emb"

task2_g2v <- emb_task2(test = test,method_name = "graph2vec",file_info = file_info,distance_type = "cosine",file_info_dups = file_info_dups)

test_ged <- fread("embeddings/ged_distance/ged_embs512_regembs_samples_all_4ep.csv")
test_ged <- test_ged[,-1]
test_ged[1] <- test_ged$emb

task2_ged <- emb_task2(test = test,method_name = "deepSNEM-GED",file_info = file_info,distance_type = "cosine",file_info_dups = file_info_dups)

test2 <- fread("embeddings/infomax/dgi_1024_tl_1.csv")
test2 <- test2[-1,-1]
colnames(test2)[1] <- "emb"

task2_IM <- emb_task2(test = test2,method_name = "deepSNEM-MI",file_info = file_info,distance_type = "cosine",file_info_dups = file_info_dups)

task2 <- bind_rows(task2_g2v,task2_IM,task2_ged)
task2$type <- factor(task2$type,levels = c("Duplicate signatures","Random signatures"))
violin_task2 <- ggplot(task2, aes(x=method, y=value, fill = type)) + 
  geom_violin(position = position_dodge(width = 1),width = 1)+geom_boxplot(position = position_dodge(width = 1),width = 0.05,
                                                                           outlier.shape = NA)+
  scale_fill_discrete(name="Embedding distance distribution",
                      labels=c("Duplicate signatures","Random signatures"))+
  ylim(0,max(task2$value))+
  xlab("")+ylab("Distance")+ 
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)

violin_task2 + theme(legend.position = "bottom")
