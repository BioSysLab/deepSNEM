library(tidyverse)
library(doFuture)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 12)
# gene baselines optimize knn on validation


evaluate_gene_perf_test_splits <- function(valsets,all,ds_path,landmark,sig_mapping,k){
  library(class)
  get_cmap_signatures <- function(cmap_path_to_gctx, sig_ids, landmark = TRUE, landmark_df = NULL) {
    
    
    library(tidyverse)
    library(cmapR)
    library(rhdf5)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    
    ds_path <- cmap_path_to_gctx
    if (landmark == TRUE) {
      
      cmap_gctx <- parse.gctx(ds_path,rid = as.character(landmark_df$`Entrez ID`), cid = sig_ids)
      cmap <- cmap_gctx@mat
      
      cmap <- cmap[as.character(landmark_df$`Entrez ID`),]
      
      rownames(cmap) <- landmark_df$Symbol
    }
    
    if (landmark == FALSE) {
      
      cmap_gctx <- parse.gctx(ds_path, cid = sig_ids)
      cmap <- cmap_gctx@mat
      
      entrez <- rownames(cmap)
      anno <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = entrez,
                                    columns = c("SYMBOL", "GENENAME","ENTREZID"),
                                    keytype = "ENTREZID")
      
      anno <- anno %>%
        filter(!is.na(SYMBOL))
      
      cmap <- cmap[anno$ENTREZID,]
      
      rownames(cmap) <- anno$SYMBOL
    }
    
    
    return(cmap)
    
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
    return(c(overall_acc,length(which(test_drugs$belongs==1))/length(drugs)))
  }
  all_sigs <- all %>% dplyr::select(sig_id,rdkit,moa_v1) %>% unique()
  all_eval <- as.data.frame(matrix(0,nrow = length(valsets),ncol=2))
  colnames(all_eval) <- c("sig_acc","drug_acc")
  for (i in 1:length(valsets)) {
    
    test_sigs <- valsets[[i]] %>% dplyr::select(sig_id,rdkit,moa_v1) %>% unique()
    train_sigs <- anti_join(all_sigs,test_sigs)
    train_sigs <- left_join(train_sigs,sig_mapping,by = c("sig_id"="sig_id2"))
    test_sigs <- left_join(test_sigs,sig_mapping,by = c("sig_id"="sig_id2"))
    genes_train <- get_cmap_signatures(ds_path,sig_ids = as.character(train_sigs$sig_id.y),landmark = T,landmark_df = landmark)
    genes_train <- t(genes_train)
    
    genes_test <- get_cmap_signatures(ds_path,sig_ids = as.character(test_sigs$sig_id.y),landmark = T,landmark_df = landmark)
    genes_test <- t(genes_test)
    set.seed(123)
    model1<- knn(train=genes_train, test=genes_test, cl=as.factor(train_sigs$moa_v1), k=k,use.all = T)
    test_sigs$predicted <- as.character(model1)
    colnames(test_sigs)[2] <- "test_rdkit"
    eval <- evaluate_predictions(predictions = test_sigs)
    all_eval[i,1] <- eval[[1]]
    all_eval[i,2] <- eval[[2]]
  }
  return(all_eval)
}
cv_knn_optimize <- function(valsets,all,ds_path,landmark,sig_mapping,k){
  results <- NULL
  results <- foreach(k_n = k) %dopar% {
    evaluate_gene_perf_test_splits(valsets, all,
                                   ds_path = ds_path,landmark = landmark,sig_mapping = sig_mapping,k = k_n)
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

# read validation and test sets
valsets <- NULL
for (i in 1:4) {
  val <- readRDS(paste0("data/gene_data/splits/rds/val_set_",i,".rds"))
  valsets[[i]] <- val 
}
test <- readRDS("data/gene_data/splits/rds/test_set.rds")

# read all sigs
all <- readRDS("data/gene_data/allsigs.rds")

# read required data files 
sig_mapping <- readRDS("data/graph_info_df/sig_mapping.rds")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "data/cmap/util_files/cmap_landmark_genes.txt")



k <- seq(1,200,5)

results <- cv_knn_optimize(valsets,all,ds_path,landmark,sig_mapping,k)

k <- 16

baseline_val <- evaluate_gene_perf_test_splits(valsets,all,
                                           ds_path = ds_path,landmark = landmark,sig_mapping = sig_mapping,k = 16)
baseline_test <- evaluate_gene_perf_test_splits(list(test),all,
                                                ds_path = ds_path,landmark = landmark,sig_mapping = sig_mapping,k = 16)

write.csv(baseline_val,"gene_baselines/knn_simple/split3/results_optimized/cv_k16.csv")
write.csv(baseline_test,"gene_baselines/knn_simple/split3/results_optimized/test_k16.csv")
