library(tidyverse)
library(doFuture)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)

evaluate_embs_test_splits <- function(valsets,all_df,all_embs,length_emb,k){
  library(class)
  evaluate_predictions <- function(predictions) {
    overall_acc <- 0
    for (i in 1:nrow(predictions)) {
      if (any(predictions$predicted[i]==predictions$moa_v1[i])) {
        overall_acc <- overall_acc+1
      }
    }
    overall_acc <- overall_acc/nrow(predictions)
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
    return(c(overall_acc,length(which(test_drugs$belongs==1))/length(drugs)))
  }
  all_sigs <- all_df %>% dplyr::select(sig_id,rdkit,moa_v1) %>% unique()
  all_eval <- as.data.frame(matrix(0,nrow = length(valsets),ncol=2))
  colnames(all_eval) <- c("sig_acc","drug_acc")
  all_embs <- left_join(all_embs,all_sigs)
  for (i in 1:length(valsets)) {
    
    test_embs <- all_embs[which(all_embs$sig_id %in% valsets[[i]]$sig_id),]
    train_embs <- anti_join(all_embs,test_embs,by = "sig_id")
    
    
    
    set.seed(123)
    model1<- knn(train=train_embs[,2:(length_emb+1)], 
                 test=test_embs[,2:(length_emb+1)], 
                 cl=as.factor(train_embs$moa_v1), k=k,use.all = T)
    #test_sigs <- valsets[[i]]
    test_embs$predicted <- as.character(model1)
    test_embs <- test_embs %>% dplyr::select(sig_id,rdkit,moa_v1,predicted)
    colnames(test_embs)[2] <- "test_rdkit"
    eval <- evaluate_predictions(predictions = test_embs)
    all_eval[i,1] <- eval[[1]]
    all_eval[i,2] <- eval[[2]]
  }
  return(all_eval)
}

cv_knn_optimize <- function(valsets,all_df,all_embs,length_emb,k){
  results <- NULL
  results <- foreach(k_n = k) %dopar% {
    evaluate_embs_test_splits(valsets = valsets,all_df=all_df,all_embs=all_embs,length_emb = length_emb,k=k_n)
  }
  k_invest <- as.data.frame(matrix(0,nrow = length(results),ncol = 2))
  for (i in 1:length(results)) {
    df <- results[[i]]
    k_invest[i,1] <- mean(df$sig_acc)
    k_invest[i,2] <- mean(df$drug_acc)
  }
  colnames(k_invest) <- c("mean_sig_acc","mean_drug_acc")
  return(k_invest)
}

splits <- readRDS("embeddings/ged_distance_semi/split3/test_splits.rds")
#val1 <- readRDS('../deepSNEM_personal/val_set_1_withoriginals.rds')
#val1 <- val1 %>% mutate(sig_id=sig_id_original)
#val2 <- readRDS('../deepSNEM_personal/val_set_2_withoriginals.rds')
#val2 <- val2 %>% mutate(sig_id=sig_id_original)
#val3 <- readRDS('../deepSNEM_personal/val_set_3_withoriginals.rds')
#val3 <- val3 %>% mutate(sig_id=sig_id_original)
#val4 <- readRDS('../deepSNEM_personal/val_set_4_withoriginals.rds')
#val4 <- val4 %>% mutate(sig_id=sig_id_original)
#valsets <- list(val1,val2,val3,val4)
valsets <- list(splits[[1]],splits[[3]],splits[[4]],splits[[5]])
test_set <- list(splits[[2]])
#test_set <- readRDS('../deepSNEM_personal/test_set_withoriginals.rds')
#test_set <- test_set %>% mutate(sig_id=sig_id_original)
all_df <- emb_proc_df
all_embs <- emb_proc

k <- seq(1,200,1)
all_df <- all_df[-which(all_df$sig_id %in% test_set$sig_id),]
all_embs <- all_embs[-which(all_embs$sig_id %in% test_set$sig_id),]
results <- cv_knn_optimize(valsets,all_df,all_embs,175,k)

k <- k[25]

baseline_val <- evaluate_embs_test_splits(valsets,all_df,all_embs,
                                               length_emb = 175,k = k)
all_df <- emb_proc_df
all_embs <- emb_proc
baseline_test <- evaluate_embs_test_splits(list(test_set),all_df,all_embs,
                                                length_emb = 175, k =k)
