library(tidyverse)
library(doFuture)
# parallel set number of workers
registerDoFuture()
plan(multiprocess,workers = 6)

# load functions
prepape_embs <- function(emb, type = "unweighted", file_info, labels){
  library(tidyverse)
  file_info <- file_info %>% dplyr::select(sig_id_original,sig_id,rdkit,cell_id,count.x)
  # remove duplicate sig id embeddings
  emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$sig_id_original)),]
  file_info <- file_info[which(as.character(file_info$sig_id_original) %in% as.character(emb$emb)),]
  # add sig id info to embs
  emb <- left_join(emb,file_info,by=c("emb"="sig_id_original"))
  # add label info to embs
  labels <- labels %>% group_by(moa_v1) %>% mutate(count = n_distinct(rdkit_graph)) %>% ungroup()
  labels <- labels %>% group_by(rdkit_graph) %>% filter(count == max(count)) %>% ungroup()
  labels <- labels %>% group_by(rdkit_graph) %>% sample_n(1) %>% ungroup()
  emb <- left_join(emb,labels,by = c("rdkit"="rdkit_graph"))
  colnames(emb)[1] <- 'sig_id_original'
  return(emb)
}
prepare_test_embs <- function(test_embs,test_df){
  test_embs <- left_join(test_embs,test_df,by = c("emb"="sig_id_original"))
  colnames(test_embs)[1] <- 'sig_id_original'
  test_embs <- test_embs %>% filter(!is.na(sig_id_original))
  return(test_embs)
}
vis_train_test_embs <- function(test_embs,train_embs,moa_n,perpl,init_dim,iter,emb_size,output_dir,name){
  library(Rtsne)
  addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    myPlot +
      guides(shape = guide_legend(override.aes = list(size = pointSize)),
             color = guide_legend(override.aes = list(size = pointSize))) +
      theme(legend.title = element_text(size = textSize), 
            legend.text  = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
  }
  
  moa <- train_embs %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
  moa_vis <- moa$moa_v1[1:moa_n]
  moa_vis <- unique(c(moa_vis,unique(as.character(test_embs$moa_v1))))
  train_embs <- train_embs[which(as.character(train_embs$moa_v1) %in% moa_vis),]
  
  train_embs$name <- ""
  test_embs <- test_embs %>% mutate(name = paste0("test_",moa_v1))
  tsne_all <- Rtsne(scale(rbind(train_embs[,2:(emb_size+1)],test_embs[,2:(emb_size+1)])), 
                    dims = 2, perplexity=perpl, 
                    verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
  
  df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2], 
                       label = as.factor(c(train_embs$moa_v1,test_embs$moa_v1)),
                       name = c(train_embs$name,test_embs$name))
  gtsne <- ggplot(df_all, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_discrete()+
    geom_text(aes(label=name),hjust=0, vjust=0,size = 0.5)
  png(file=paste0(output_dir,"/",name,".png"),width=9,height=9,units = "in",res=600)
  print(addSmallLegend(gtsne))
  dev.off()
  
  
}
predict_test_moa_emb <- function(test_embs,train_embs,emb_size,k){
  library(tidyverse)
  library(class)
  train_embs <- train_embs %>% filter(!is.na(moa_v1))
  set.seed(123)
  model1<- knn(train=train_embs[,2:(emb_size+1)], test=test_embs[,2:(emb_size+1)], cl=as.factor(train_embs$moa_v1), k=k,use.all = T)
  test_embs$predicted <- as.character(model1)
  results <- test_embs %>% dplyr::select(test_rdkit,sig_id,sig_id_original,moa_v1,predicted,cell_id,pert_iname) %>% unique()
  return(results)
}
knn_dist_predictions <- function(train_embs,test_embs,train_labels,k,emb_length){
  cosine_dist <- function(a,b){
    a <- as.vector(as.matrix(a))
    b <- as.vector(as.matrix(b))
    cos_dist <- 1 - (a%*%b)/(sqrt(a%*%a)*sqrt(b%*%b))
    return(cos_dist)
  }
  get_k_nearest <- function(a,k,col_names,train_embs){
    a <- as.vector(as.matrix(a))
    names(a) <- col_names
    a <- sort(a,decreasing = F)
    neighbors <- names(a)[1:k]
    moas <- as.character(train_embs$moa_v1[which(train_embs$sig_id_original %in% neighbors)])
    un_moas <- unique(moas)
    freq <- NULL
    for (j in 1:length(un_moas)) {
      freq[j] <- length(which(moas %in% un_moas[j]))/length(moas)
    }
    names(freq) <- un_moas
    assigned_moa <- names(freq)[which(freq==max(freq))]
    return(assigned_moa)
  }
  d1 <- apply(test_embs[,2:(emb_length+1)], 1, function(x) apply(train_embs[,2:(emb_length+1)], 1, function(y) cosine_dist(x,y)))
  d1 <- t(d1)
  rownames(d1) <- as.character(test_embs$sig_id_original)
  colnames(d1) <- as.character(train_embs$sig_id_original)
  moas <- apply(d1,1,get_k_nearest,k,colnames(d1),train_embs)
  results <- test_embs %>% dplyr::select(test_rdkit,sig_id,sig_id_original,moa_v1,cell_id,pert_iname) %>% unique()
  results$predicted <- moas
  return(results)
}
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
  return(list(overall_acc,test_drugs,length(which(test_drugs$belongs==1))/length(drugs),assigned_moa))
}
investigate_k <- function(k,test_embs,train_embs,emb_length,knn_type){
  results <- matrix(666,nrow = length(k),ncol = 3)
  for (i in 1:length(k)) {
    if (knn_type == "normal") {
      predictions <- predict_test_moa_emb(test_embs = test_embs,train_embs = train_embs,
                                          emb_size = emb_length,
                                          k = k[i])
    }
    if (knn_type == "modified") {
      predictions <- knn_dist_predictions(train_embs = train_embs,test_embs = test_embs,
                                          train_labels = train_embs$moa_v1,k = k[i],emb_length = emb_length)
    }
    eval <- evaluate_predictions(predictions = predictions)
    results[i,1] <- k[i]
    results[i,2] <- eval[[1]]
    results[i,3] <- eval[[3]]
  }
  results <- as.data.frame(results)
  colnames(results) <- c("k","sig_acc","drug_acc")
  return(results)
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

#####
# load required files for the functions

file_info <- readRDS("data/graph_info_df/file_info_nodups.rds") %>% select(-emb,-files_combined) %>% unique()
file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds") %>% select(-emb,-files_combined,-sig_id.y) %>% unique()
##na valo mesa sig_id_original apo to file_info_weighted
info_origins <- readRDS("data/graph_info_df/file_info_weighted.rds") %>% select(sig_id,sig_id_original) %>% unique()
file_info <- left_join(file_info,info_origins,by="sig_id")
file_info_dups <- left_join(file_info_dups,info_origins,by="sig_id")
labels <- readRDS("data/cmap/labels/labels_first_pass.rds")
allpairs <- readRDS("embeddings/ged_distance_semi/split3/allpairs3_genes.rds")
sigs <- unique(c(as.character(allpairs$sig_id_original.x),as.character(allpairs$sig_id_original.y))) # all graphs with moa

# load train (ALL), test, validation  embeddings test, validation dataframes

emb <- read.csv("../deepSNEM_personal/tf_embs_q1.csv") # ALL
test_embs <- read.csv("../deepSNEM_personal/Ugraphemp/embeddings/gsea_embs512_best_genes_test_0.csv") # TEST
val_embs <- read.csv("../deepSNEM_personal/Ugraphemp/embeddings/gsea_embs512_best_genes_val_0.csv") # VAL

# load test,val dataframe 

test_df <- readRDS("embeddings/ged_distance_semi/split3/test_set_withoriginals.rds")
test_df <- test_df %>% select(-graphs,-emb) %>% unique()
val_df <- readRDS("embeddings/ged_distance_semi/split3/val_set_1_withoriginals.rds")
val_df <- val_df %>% select(-graphs,-emb) %>% unique()

# preprocess embeddings, for each sig id keep one embedding (either averaged or random selection)
emb <- emb[,-1]
colnames(emb)[1] <- "emb"
emb_proc <- prepape_embs(emb = emb,file_info = file_info,labels = labels)

# preprocess test and val dataframes (keep only graphs that were present in the semisupervised task)

colnames(test_df) <- c("test_rdkit","sig_id","pert_iname","cell_id","moa_v1","sig_id_original")
test_df <- test_df[which(test_df$sig_id_original %in% sigs),]
colnames(val_df) <- c("test_rdkit","sig_id","pert_iname","cell_id","moa_v1","sig_id_original")
val_df <- val_df[which(val_df$sig_id_original %in% sigs),]

# prepare test and validation embeddings, one embedding for each sig id

test_embs <- test_embs[,-1]
colnames(test_embs)[1] <- "emb"
val_embs <- val_embs[,-1]
colnames(val_embs)[1] <- "emb"
test_embs <- prepare_test_embs(test_embs = test_embs,test_df = test_df)
val_embs <- prepare_test_embs(test_embs = val_embs,test_df = val_df)

# remove from the train embs the test and val embs to run KNN
train_embs <- emb_proc[-which(emb_proc$sig_id_original %in% test_embs$sig_id_original),]
train_embs <- train_embs[-which(train_embs$sig_id_original %in% val_embs$sig_id_original),]
# remove the train embeddings without moa 
train_embs <- train_embs %>% filter(moa_v1 != "") %>% filter(!is.na(moa_v1))


#visualize test/val and train embeddings
vis_train_test_embs(test_embs = test_embs,train_embs = train_embs,moa_n = 10,perpl = 5,init_dim = 80,iter = 1000,emb_size = 256,
                    output_dir = "vis",name = "test_set_emb512_genes_2ep"  )

# investigate k on the validation embeddings on KNN with cosine or euclidian distance

k <- seq(1,200,5)

results_modified <- NULL
results_modified <- foreach(k_n = k) %dopar% {
  investigate_k(k = k_n,test_embs = test_embs,train_embs = train_embs,emb_length = 512,knn_type = "modified")
}
results_df <- results_modified[[1]]
for (i in 2:length(results_modified)) {
  results_df <- rbind(results_df,results_modified[[i]])
}

# investigate k on KNN with euclidian distance
#results_normal <- NULL
#results_normal <- foreach(k_n = k) %dopar% {
#  investigate_k(k = k_n,test_embs = test_embs,train_embs = train_embs,emb_length = 512,knn_type = "normal")
#}
#results_df_normal <- results_normal[[1]]
#for (i in 2:length(results_normal)) {
#  results_df_normal <- rbind(results_df_normal,results_normal[[i]])
#}

# KNN with test predictions with k from the validation set
k <- 61
# predictions from cosine KNN
predictions <- knn_dist_predictions(train_embs = train_embs,test_embs = test_embs,train_labels = train_embs$moa_v1,
                                    emb_length = 512,k = k)

# predictions from euclidian KNN
#predictions <- predict_test_moa_emb(test_embs = test_embs,train_embs = train_embs,emb_size = 512,k=k)
# predictions from euclidian KNN
# evaluate predictions
eval_emb <- evaluate_predictions(predictions = predictions)

write_results(eval_emb,output_dir = "embeddings/ged_distance_semi/split3/results",test_name = "genes_semi_emb512_test_2ep")




