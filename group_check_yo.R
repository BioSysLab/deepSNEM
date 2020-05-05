group_check <- function(group,group_name,original_big_group,
                        labels,emb_test,file_info,emb_size,output_dir,
                        label_space_only, tsne_perpl, umap_n,
                        cell_specific,cell_line,genes, ds_path, 
                        landmark,sig_map,tsne_perpl_genes, init_dim, init_dim_genes){
  library(tidyverse)
  library(clustree)
  library(Rtsne)
  library(ggplot2)
  library(umap)
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
  dir.create(output_dir,recursive = T)
  labels$moa_v1 <- labels$moa
  labels$moa_v1[which(labels$moa %in% group)] <- group_name
  labels_v1 <- left_join(labels,file_info, by = c("rdkit_graph"="rdkit"))
  
  labels_v1 <- labels_v1 %>% mutate(test_labels = if_else(condition = (moa_v1 == group_name),true = group_name,false = "nada"))
  
  if (original_big_group != "") {
    labels_v1 <- labels_v1 %>% 
      mutate(test_labels = if_else(condition = (moa == original_big_group),true = original_big_group,false = test_labels))
  }
  if (label_space_only) {
    emb_test <- emb[which(emb$emb %in% labels_v1$emb),]
  }
  emb_test <- left_join(emb_test,labels_v1,by="emb")
  if (cell_specific) {
    emb_test <- emb_test[which(as.character(emb_test$cell_id) %in% cell_line),]
  }
  print(nrow(emb_test))
  tsne_test <- Rtsne(scale(emb_test[,2:(emb_size+1)]), dims = 2, perplexity=tsne_perpl, verbose=TRUE, max_iter = 1000,initial_dims = init_dim)
  df_test <- data.frame(V1 = tsne_test$Y[,1], V2 =tsne_test$Y[,2], label = as.factor(emb_test$test_labels))
  colors <- c( "#E69F00","#999999", "#56B4E9")
  names(colors) <- c(as.character(original_big_group),"nada",as.character(group_name))
  gtsne <- ggplot(df_test, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_manual(values=colors)
  png(file=paste0(output_dir,"/",as.character(group_name),"tsne.png"),width=9,height=9,units = "in",res=300)
  print(gtsne)
  dev.off()
  umap.defaults$n_neighbors <- umap_n
  map_test <- umap(scale(emb_test[,2:(emb_size+1)]))
  df_map <- data.frame(V1 = map_test$layout[,1], V2 = map_test$layout[,2], label = as.factor(emb_test$test_labels))
  gmap <- ggplot(df_map, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_manual(values=colors)
  png(file=paste0(output_dir,"/",as.character(group_name),"umap.png"),width=9,height=9,units = "in",res=300)
  print(gmap)
  dev.off()
  emb_test <- left_join(emb_test,file_info,by = "emb")
  emb_test <- left_join(emb_test,sig_map,by=c("sig_id.y"="sig_id2"))
  if (genes) {
    emb_genes <- get_cmap_signatures(as.character(emb_test$sig_id),cmap_path_to_gctx = ds_path,landmark = T, landmark_df = landmark)
    emb_genes <- t(emb_genes)
    emb_genes <- as.data.frame(emb_genes) %>% rownames_to_column("sig_id")
    tsne_genes <- Rtsne(scale(emb_genes[,2:ncol(emb_genes)]), dims = 2, perplexity=tsne_perpl_genes, verbose=TRUE, max_iter = 1000,initial_dims = init_dim_genes)
    df_tsne_genes <- data.frame(V1 = tsne_genes$Y[,1], V2 =tsne_genes$Y[,2], label = as.factor(emb_test$test_labels))
    gtsne_genes <- ggplot(df_tsne_genes, aes(V1, V2))+
      geom_point(aes(color = label),show.legend = T) + scale_color_manual(values=colors)
    #map_genes <- umap(scale(emb_genes[,2:ncol(emb_genes)]))
    #df_map_genes <- data.frame(V1 = map_genes$layout[,1], V2 = map_genes$layout[,2], label = as.factor(emb_test$test_labels))
    #gmap_genes <- ggplot(df_map_genes, aes(V1, V2))+
      #geom_point(aes(color = label),show.legend = T) + scale_color_manual(values=colors)
    png(file=paste0(output_dir,"/",as.character(group_name),"tsne_genes.png"),width=9,height=9,units = "in",res=300)
    print(gtsne_genes)
    dev.off()
    #png(file=paste0(output_dir,"/",as.character(group_name),"umap_genes.png"),width=9,height=9,units = "in",res=300)
    #print(gmap_genes)
    #dev.off()
  }
}

# file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)
# embs
emb <- read.csv("embeddings/graph2vec/emb_activity_1_epoch.csv")
colnames(emb)[1] <- "emb"
#emb <- emb %>% mutate(emb = str_remove_all(string = emb,pattern = ".csv"))
# keep only 1 emb per sig id for clustering
# keep the respective embeddings

emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
file_info <- file_info[which(as.character(file_info$emb) %in% as.character(emb$emb)),]

file_info <- file_info %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]


labels <- readRDS("data/cmap/labels.rds")

dna_rna_damage <- c("topoisomerase inhibitor","RNA synthesis inhibitor|topoisomerase inhibitor",
                    "chelating agent|topoisomerase inhibitor","RNA synthesis inhibitor","DNA synthesis inhibitor",
                    "DNA alkylating agent","DNA alkylating agent|DNA synthesis inhibitor","DNA inhibitor",
                    "DNA replication inhibitor|STAT inhibitor")

group_name <- "dna_rna_damage"

output_dir <- paste0("group_check/",as.character(group_name))

obg = "topoisomerase inhibitor"

group_check(group = dna_rna_damage,group_name = "dna_rna_damage",original_big_group = obg,
            labels = labels, emb_test = emb,file_info = file_info,emb_size = 128,
            output_dir = output_dir, label_space_only = F, tsne_perpl = 5, init_dim = 50,umap_n = 50, cell_specific = T,cell_line = "A375",
            genes = T,
            ds_path = ds_path, landmark = landmark, sig_map = sig_map, tsne_perpl_genes = 5, init_dim_genes = 50)
