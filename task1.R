emb_task1 <- function(test, method_name, file_info, distance_type){
  # test is a dataframe where there is a column emb which is in the
  # form of the column emb of file info
  # file info, dataframe 
  # distance, character string, either "cosine" or "euclidian"
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
  
  # Task 1 same sig id similar graph embeddings
  
  file_info_t1 <- file_info %>% filter(count_per_sig > 1)
  sigs <- unique(as.character(file_info_t1$sig_id))
  
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
  
  sampled_sigs <- sample(sigs,100)
  all_sampled_dists <- NULL
  file_info_t1$sig_id <- as.character(file_info_t1$sig_id)
  for (i in 1:length(sampled_sigs)) {
    filt <- file_info_t1 %>% filter(sig_id == sampled_sigs[i])
    df <- test[which(test$emb %in% filt$emb),]
    sig_dists <- distance_function(df = df)
    all_sampled_dists <- rbind(sig_dists,all_sampled_dists)
  }
  
  sigs2 <- unique(as.character(file_info$sig_id))
  set.seed(1337)
  sigs_random <- sample(sigs2,90)
  
  random_embs <- NULL
  for (i in 1:length(sigs_random)) {
    filt <- file_info %>% filter(sig_id == sigs_random[i])
    random_emb <- sample_n(test[which(test$emb %in% filt$emb),],1)
    random_embs <- rbind(random_emb,random_embs)
  }
  all_random_dists <- distance_function(random_embs)
  
  all_sampled_dists$method <- method_name
  all_sampled_dists$type <- "Same signatures"
  
  all_random_dists$method <- method_name
  all_random_dists$type <- "Different signatures"
  
  task1 <- bind_rows(all_sampled_dists,all_random_dists)
  task1$value <- task1$value/max(task1$value)
  return(task1)
}

# file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
library(data.table)
test <- fread("embeddings/graph2vec/emb_activity_1_epoch.csv")
colnames(test)[1] <- "emb"
task1_g2v <- emb_task1(test = test,method_name = "graph2vec",file_info = file_info,distance_type = "cosine")

test2 <- fread("embeddings/infomax/dgi_1024_tl_1.csv")
test2 <- test2[-1,-1]
colnames(test2)[1] <- "emb"

task1_IM <- emb_task1(test = test2,method_name = "deepSNEM-MI",file_info = file_info,distance_type = "cosine")

test_ged <- fread("embeddings/ged_distance/ged_embs512_regembs_samples_all_4ep.csv")
test_ged <- test_ged[,-1]
test_ged[1] <- test_ged$emb

task1_ged <- emb_task1(test = test_ged, method_name = "deepSNEM-GED", file_info = file_info, distance_type = "cosine")

task1 <- bind_rows(task1_g2v,task1_IM,task1_ged)
task1$type <- factor(task1$type,levels = c("Same signatures","Different signatures"))
violin_task1 <- ggplot(task1, aes(x=method, y=value, fill = type)) + 
  geom_violin(position = position_dodge(width = 1),width = 1)+geom_boxplot(position = position_dodge(width = 1),width = 0.05,
                                                                             outlier.shape = NA)+
  scale_fill_discrete(name="Embedding distance distribution",
                      labels=c("Same signatures","Different signatures"))+
  ylim(0,max(task1$value))+
  xlab("")+ylab("Distance")+ 
  theme(axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(family = "serif",size = 20),legend.position = "bottom")+
  theme_minimal(base_family = "serif",base_size = 20)

violin_task1 + theme(legend.position = "bottom")

