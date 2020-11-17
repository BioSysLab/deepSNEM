library(tidyverse)

allpairs <- readRDS("data/graph_additional/pairs/splits/split3/allpairs3.rds")
length(which(allpairs$sig_id.x==allpairs$sig_id.y))
length(which(allpairs$label==1))
allpairs <- allpairs %>% filter(sig_id.x != sig_id.y)
allpairs <- allpairs %>% filter(moa_v1.x != "") %>% filter(moa_v1.y != "") %>% filter(!is.na(moa_v1.x)) %>% filter(!is.na(moa_v1.y))
length(unique(c(allpairs$graph.x,allpairs$graph.y)))
length(unique(c(allpairs$sig_id.x,allpairs$sig_id.y)))
test <- allpairs %>% group_by(sig_id.x,sig_id.y) %>% summarise(count = n())

# read file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")

graphs <- readRDS("data/graph_additional/pairs/graphs_no_dupl.rds")

# turn to sig id instead of graphs

sigs <- graphs %>% dplyr::select(sig_id,pert_iname,cell_id,rdkit) %>% unique()

# add labels to sigs

labels <- readRDS("data/cmap/labels/labels_first_pass.rds")

sigs <- left_join(sigs,labels,by=c("rdkit"="rdkit_graph"))
# remove the sigs with no label
sigs <- sigs %>% filter(!is.na(moa_v1)) %>% filter(moa_v1 != "")
# remove from the splitting process the drugs with multiple moas
sigs <- sigs %>% group_by(rdkit) %>% mutate(duplicate_moa = n_distinct(moa_v1)) %>% ungroup()
sigs <- sigs %>% filter(duplicate_moa==1)
# count up unique sig ids per rdkit

sigs <- sigs %>% group_by(rdkit) %>% mutate(rdkit_count = n_distinct(sig_id)) %>% ungroup()

# count up number of unique drugs for each moa

sigs <- sigs %>% group_by(moa_v1) %>% mutate(moa_count = n_distinct(sig_id)) %>% ungroup()
sigs <- sigs %>% group_by(moa_v1) %>% mutate(moa_count2 = n_distinct(rdkit)) %>% ungroup()

moa <- sigs %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))




min_moa <- 25
n_split <- 5
get_test_drugs <- function(n_split,sigs,min_moa){
  candidates <- sigs %>% dplyr::select(rdkit,rdkit_count,moa_v1,moa_count,moa_count2) %>% unique()
  candidates <- candidates %>% filter(moa_count>=min_moa) %>% filter(moa_count2 >= n_split)
  test_list <- list(0)
  sigs <- sigs %>% dplyr::select(sig_id,cell_id,pert_iname,rdkit)
  for (i in 1:n_split) {
    random_test <- candidates %>% group_by(moa_v1) %>% sample_n(1) %>% ungroup()
    candidates <- candidates[-which(candidates$rdkit%in%random_test$rdkit),]
    random_test <- left_join(random_test,sigs,by=c("rdkit"="rdkit"))
    test_list[[i]] <- random_test
  }
  return(test_list)
}


test_splits <- get_test_drugs(5,sigs,25)
test_splits <- readRDS("data/graph_additional/pairs/splits/split3/test_splits.rds")
saveRDS(test_splits,"data/graph_additional/pairs/splits/split3/test_splits.rds")
evaluate_gene_perf_test_splits <- function(test_splits,sigs,ds_path,landmark,sig_mapping,k){
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
    return(list(overall_acc,test_drugs,length(which(test_drugs$belongs==1))/length(drugs),assigned_moa))
  }
  all_sigs <- sigs %>% dplyr::select(sig_id,rdkit,moa_v1) %>% unique()
  all_eval <- list(0)
  for (i in 1:length(test_splits)) {
    
    test_sigs <- test_splits[[i]] %>% dplyr::select(sig_id,rdkit,moa_v1) %>% unique()
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
    all_eval[[i]] <- eval
  }
  return(all_eval)
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

sig_mapping <- readRDS("data/graph_info_df/sig_mapping.rds")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "data/cmap/util_files/cmap_landmark_genes.txt")

baseline <- evaluate_gene_perf_test_splits(test_splits = test_splits,sigs = sigs,
                                           ds_path = ds_path,landmark = landmark,sig_mapping = sig_mapping,k = 3)

write_results(eval = baseline[[1]],output_dir = "gene_baselines/knn_simple/split3/results",test_name = "val_set_1")
write_results(eval = baseline[[2]],output_dir = "gene_baselines/knn_simple/split3/results",test_name = "test_set")

test_sigs <- test_splits[[2]]
test_embs <- left_join(test_sigs,graphs,by=c("sig_id"="sig_id"))
test_embs <- test_embs %>% dplyr::select(rdkit.x,sig_id,graphs,pert_iname.x,cell_id.x,emb,moa_v1)
saveRDS(test_embs,"data/graph_additional/pairs/splits/split3/test_set.rds")
write.csv(test_embs,"data/graph_additional/pairs/splits/split3/test_set.csv")
test_set_number <- 2
for (i in 1:(length(test_splits))) {
  if (i != test_set_number) {
    test_sigs <- test_splits[[1]]
    test_embs <- left_join(test_sigs,graphs,by=c("sig_id"="sig_id"))
    test_embs <- test_embs %>% dplyr::select(rdkit.x,sig_id,graphs,pert_iname.x,cell_id.x,emb,moa_v1)
    saveRDS(test_embs,paste0("data/graph_additional/pairs/splits/split3/val_set_",i,".rds"))
    write.csv(test_embs,paste0("data/graph_additional/pairs/splits//split3/val_set_",i,".csv"))
  }
}

# reset all pairs labels

allpairs <- allpairs %>% dplyr::select(-label,-label1,-label2,-label3,-label4)
allpairs$label <- 0
allpairs$label[allpairs$moa_v1.x==allpairs$moa_v1.y] <- 1
length(which(allpairs$label==1))
#allpairs$label[allpairs$identifier.x==allpairs$identifier.y] <- 1
length(which(allpairs$label==1))
test_id <- unique(c(which(as.character(allpairs$sig_id.x) %in% test_sigs$sig_id),
                    which(as.character(allpairs$sig_id.y) %in% test_sigs$sig_id)))

#set labels of test all to 0
allpairs$label[test_id] <- 0

allpairs$label1 <- allpairs$label
allpairs$label2 <- allpairs$label
allpairs$label3 <- allpairs$label
allpairs$label4 <- allpairs$label
#allpairs$label5 <- allpairs$label

val_1 <- test_splits[[1]]
val_2 <- test_splits[[2]]
val_3 <- test_splits[[3]]
val_4 <- test_splits[[5]]
#val_5 <- test_splits[[5]]

val_id_1 <- unique(c(which(as.character(allpairs$sig_id.x) %in% val_1$sig_id),
                     which(as.character(allpairs$sig_id.y) %in% val_1$sig_id)))
val_id_2 <- unique(c(which(allpairs$rdkit.x %in% val_2$rdkit),which(allpairs$rdkit.y %in% val_2$rdkit)))
val_id_3 <- unique(c(which(allpairs$rdkit.x %in% val_3$rdkit),which(allpairs$rdkit.y %in% val_3$rdkit)))
val_id_4 <- unique(c(which(allpairs$rdkit.x %in% val_4$rdkit),which(allpairs$rdkit.y %in% val_4$rdkit)))
#val_id_5 <- unique(c(which(allpairs$rdkit.x %in% val_5$rdkit),which(allpairs$rdkit.y %in% val_5$rdkit)))

allpairs$label1[val_id_1] <- 0
allpairs$label2[val_id_2] <- 0
allpairs$label3[val_id_3] <- 0
allpairs$label4[val_id_4] <- 0
#allpairs$label5[val_id_5] <- 0

saveRDS(allpairs3,"data/graph_additional/pairs/splits/split3/allpairs3.rds")
write.csv(allpairs3,"data/graph_additional/pairs/splits/split3/allpairs3.csv")

#saveRDS(test_splits,"data/graph_additional/pairs/splits/split2/test_splits.rds")

labeled <- allpairs %>% filter(label1 == 1)
length(which(labeled$sig_id.x %in% test_splits[[1]]$sig_id))

pairs_rem <- anti_join(allpairs2,allpairs)
pairs_rem <- pairs_rem %>% dplyr::select(-label,-label1,-label2,-label3,-label4)
pairs_rem$label <- 0
pairs_rem$label1 <- 0 

pairs_rem2 <- pairs_rem %>% filter(!is.na(identifier.x) & !is.na(identifier.y))
pairs_rem <- anti_join(pairs_rem,pairs_rem2)
pairs_rem2$label[pairs_rem2$identifier.x==pairs_rem2$identifier.y] <- 1
pairs_rem2$label1 <- pairs_rem2$label

labeled_rem2 <- pairs_rem2 %>% filter(label==1)

pairs_rem <- rbind(pairs_rem,pairs_rem2)

allpairs3 <- rbind(allpairs3,pairs_rem)
