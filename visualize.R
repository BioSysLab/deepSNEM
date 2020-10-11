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

file_info <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/graph_info_df/file_info_nodups.rds")
file_info_dups <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/graph_info_df/file_info_dups.rds")
emb <- read.csv("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/embeddings/deep_graph_infomax/semi_supervised/DGI_JSD_512_GO_semi.csv")
#test_files <- as.character(emb$X)
#test_files <- as.data.frame(test_files)
#test_files <- left_join(test_files,file_info,by = c("test_files"="files_combined"))
#test_files <- test_files %>% dplyr::select(test_files,emb)
#test_files <- left_join(test_files,file_info_dups,by = c("test_files"="files_combined"))
#test_files <- test_files %>% dplyr::select(test_files,emb.x,emb.y) %>% mutate(emb = if_else(condition = is.na(emb.x),true = emb.y,false = emb.x))
#test_files <- test_files %>% dplyr::select(test_files,emb)
#emb <- left_join(test_files,emb,by= c("test_files"="X"))
emb <- emb[,-1]
colnames(emb)[1] <- "emb"
#emb <- readRDS("embeddings/ged_distance/ged_embs512_seen_4ep.rds")
labels <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/cmap/labels/labels_first_pass.rds")
emb_proc <- prepape_embs(emb = emb,file_info = file_info,labels = labels,keep_one = F ,ave = T,n_emb = (ncol(emb)-1))
#emb <- emb_proc
visualize_moa_emb <- function(emb,output_dir,moa_n,emb_size,perpl_emb,iter,init_dim,name,scale = T){
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
  emb <- emb %>% filter(!is.na(moa_v1))
  moa <- emb %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
  moa_vis <- moa$moa_v1[1:moa_n]
  emb <- emb[which(as.character(emb$moa_v1) %in% moa_vis),]
  if (scale) {
    tsne_emb <- Rtsne(scale(emb[,2:(emb_size+1)]), dims = 2, perplexity=perpl_emb, verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
  } else {
    tsne_emb <- Rtsne(emb[,2:(emb_size+1)], dims = 2, perplexity=perpl_emb, verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
  }
  df_emb <- data.frame(V1 = tsne_emb$Y[,1], V2 =tsne_emb$Y[,2], label = as.factor(emb$moa_v1))
  gtsne <- ggplot(df_emb, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_discrete()
  png(file=paste0(output_dir,"/",as.character(name),"_moa_tsne.png"),width=9,height=9,units = "in",res=300)
  print(addSmallLegend(gtsne))
  dev.off()
}

visualize_moa_emb(emb_proc,output_dir = "/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/vis",moa_n = 10,emb_size = 1024,
              perpl_emb = 5,init_dim = 80,iter = 2000,name = "vis_dgi_jsd_512_semi", scale = F)


visualize_moa_genes <- function(emb,output_dir,moa_n,perpl,iter,init_dim,name,ds_path,landmark,sig_map){
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
  moa <- emb %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
  moa_vis <- moa$moa_v1[1:moa_n]
  emb <- emb[which(as.character(emb$moa_v1) %in% moa_vis),]
  emb <- left_join(emb,sig_map,by = c("sig_id"="sig_id2"))
  genes <- get_cmap_signatures(ds_path,sig_ids = as.character(emb$sig_id.y),landmark = T,landmark_df = landmark)
  genes <- t(genes)
  tsne_emb <- Rtsne(scale(genes), dims = 2, perplexity=perpl, verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
  df_emb <- data.frame(V1 = tsne_emb$Y[,1], V2 =tsne_emb$Y[,2], label = as.factor(emb$moa_v1))
  gtsne <- ggplot(df_emb, aes(V1, V2))+
    geom_point(aes(color = label),show.legend = T) + scale_color_discrete()
  png(file=paste0(output_dir,"/",as.character(name),"_moa_tsne.png"),width=9,height=9,units = "in",res=300)
  print(addSmallLegend(gtsne))
  dev.off()
  return(genes)
}
sig_mapping <- readRDS("data/graph_info_df/sig_mapping.rds")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "data/cmap/util_files/cmap_landmark_genes.txt")

visualize_moa_genes(emb_proc,output_dir = "vis",moa_n = 10,
                  perpl = 15,init_dim = 80,iter = 1000,name = "test_transformer_256_tl1_leaky_relu_mean",
                  ds_path = NULL, landmark = NULL,sig_map = sig_mapping)




