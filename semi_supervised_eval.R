library(tidyverse)

prepape_embs <- function(emb, type = "unweighted", file_info, labels, keep_one, ave, n_emb){
  library(tidyverse)
  file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)
  # remove duplicate sig id embeddings
  emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
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
    file_info_1 <- file_info %>% dplyr::select(sig_id,rdkit) %>% unique()
    emb <- left_join(aver,file_info_1,by = c("Group.1"="sig_id"))
    emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
    colnames(emb)[1] <- "sig_id"
  }
  return(emb)
}

prepare_test_embs <- function(test_embs,test_df,ave_drug,ave_sig,emb_n){
  test_embs <- left_join(test_embs,test_df,by = "emb")
  if (ave_sig) {
    aver <- aggregate(test_embs[, 2:(emb_n+1)], list(test_embs$sig_id), mean)
    test_df2 <- test_df %>% dplyr::select(sig_id,test_rdkit,cell_id,moa_v1,pert_iname) %>% unique()
    test_embs <- left_join(aver,test_df2,by = c("Group.1"="sig_id"))
    colnames(test_embs)[1] <- "sig_id"
  }
  return(test_embs)
}

predict_test_moa_emb <- function(test_embs,train_embs,moa_n,perpl,init_dim,iter,emb_size,output_dir,vis = F,k){
  library(tidyverse)
  library(Rtsne)
  library(class)
  addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_text(size = textSize), 
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }
  cos_dist <- function(a,b) {
    return(1-sum(a*b)/sqrt(sum(a^2)*sum(b^2)) )
  }
  
  train_embs <- train_embs %>% filter(!is.na(moa_v1))
  moa <- train_embs %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
  moa_vis <- moa$moa_v1[1:moa_n]
  train_embs <- train_embs[which(as.character(train_embs$moa_v1) %in% moa_vis),]
  
  train_embs$name <- ""
  test_embs <- test_embs %>% mutate(name = paste0("test_",moa_v1))
  
  if (vis) {
    tsne_all <- Rtsne(scale(rbind(train_embs[,2:(emb_size+1)],test_embs[,2:(emb_size+1)])), 
                      dims = 2, perplexity=perpl, 
                      verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
    
    df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2], 
                         label = as.factor(c(train_embs$moa_v1,test_embs$moa_v1)),
                         name = c(train_embs$name,test_embs$name))
    gtsne <- ggplot(df_all, aes(V1, V2))+
      geom_point(aes(color = label),show.legend = T) + scale_color_discrete()+
      geom_text(aes(label=name),hjust=0, vjust=0,size = 0.5)
    png(file=paste0(output_dir,"/test_set_1.png"),width=9,height=9,units = "in",res=600)
    print(addSmallLegend(gtsne))
    dev.off()
  }
  set.seed(123)
  model1<- knn(train=train_embs[,2:513], test=test_embs[,2:513], cl=as.factor(train_embs$moa_v1), k=k,use.all = T)
  test_embs$predicted <- as.character(model1)
  results <- test_embs %>% dplyr::select(test_rdkit,sig_id,moa_v1,predicted,cell_id,pert_iname) %>% unique()
  return(results)
}

evaluate_predictions <- function(predictions) {
  overall_acc <- length(which(predictions$predicted==predictions$moa_v1))/nrow(predictions)
  drugs <- unique(as.character(predictions$test_rdkit))
  assigned_moa <- list(0)
  for (i in 1:length(drugs)) {
    filt <- predictions[which(predictions$test_rdkit %in% drugs[i]),]
    moas <- as.character(unique(filt$predicted))
    freq <- NULL
    for (j in 1:length(moas)) {
      freq[j] <- length(which(filt$predicted %in% moas[j]))/nrow(filt)
    }
    names(freq) <- moas
    assigned_moa_n <- names(freq)[which(freq==max(freq))]
    assigned_moa[[i]] <- assigned_moa_n
  }
  
  test_drugs <- as.data.frame(drugs)
  test_drugs <- left_join(test_drugs,predictions,by = c("drugs"="test_rdkit")) %>% dplyr::select(drugs,moa_v1) %>% unique()
  test_drugs$belongs <- 0
  for (i in 1:nrow(test_drugs)) {
    if (any(test_drugs$moa_v1[i]==assigned_moa[[i]])) {
      test_drugs$belongs[i] <- 1
    }
  }
  return(list(overall_acc,test_drugs,length(which(test_drugs$belongs==1))/length(drugs),assigned_moa))
}

write_results <- function(eval,output_dir,test_name){
  sig_accuracy <- eval[[1]]
  drug_accuracy <- eval[[3]]
  df <- eval[[2]]
  assigned <- eval[[4]]
  df$assigned <- assigned
  df$sig_accuracy <- sig_accuracy
  df$drug_accuracy <- drug_accuracy
  df <- df %>% unnest(assigned)
  df <- df %>% group_by(moa_v1) %>% mutate(assigned = paste(assigned,collapse = "/")) %>% ungroup() %>% unique()
  write.csv(df,paste0(output_dir,"/",test_name,"_eval_results.csv"))
}



file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
emb <- read.csv("embeddings/ged_distance_semi/split2/ss_1ep_512_cosine_train_split2.csv")
emb <- emb[,-1]
colnames(emb)[1] <- "emb"
labels <- readRDS("data/cmap/labels/labels_first_pass.rds")
emb_proc <- prepape_embs(emb = emb,file_info = file_info,labels = labels,keep_one = F ,ave = T,n_emb = (ncol(emb)-1))



test_df <- readRDS("data/graph_additional/pairs/splits/split2/test_set.rds")
colnames(test_df) <- c("test_rdkit","sig_id","files_combined","pert_iname","cell_id","emb","moa_v1")
test_embs <- read.csv("embeddings/ged_distance_semi/split2/1epoch/ged_ss_1ep_test_set_split2.csv")
test_embs <- test_embs[,-1]
colnames(test_embs)[1] <- "emb"

train_embs <- emb_proc
test_embs <- prepare_test_embs(test_embs = test_embs,test_df = test_df,ave_drug = F,ave_sig = T,emb_n=512)

predictions <- predict_test_moa_emb(test_embs = test_embs,train_embs = train_embs,
                                moa_n = 10,perpl = 5,init_dim = 80,iter = 2000,
                                emb_size = 512,output_dir = "",vis = F,
                                k = 1)


eval_emb <- evaluate_predictions(predictions = predictions)

write_results(eval_emb,output_dir = "embeddings/ged_distance_semi/split2/1epoch/results",test_name = "test_set")


