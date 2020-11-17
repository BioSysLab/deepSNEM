eval_emb <- function(test, file_info, distance_type, file_info_dups, ds_path,landmark,sig_mapping, sims, output_dir, task3, genes){
  # test is a dataframe where there is a column emb which is in the
  # form of the column emb of file info
  # file info, dataframe 
  # distance, character string, either "cosine" or "euclidian"
  library(tidyverse)
  library(lsa)
  library(Rtsne)
  library(umap)
  # keep some columns of file_info
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
  random_sigs <- sample(sigs,round(sqrt(2*nrow(all_sampled_dists))))
  all_random_dists <- NULL
  ids_random_sigs <- which(file_info_t1$sig_id %in% random_sigs)
  random_info <- file_info_t1[ids_random_sigs,]
  random_info <- random_info %>% group_by(sig_id) %>% filter(emb == emb[1])
  ids_random_emb <- which(test$emb %in% random_info$emb)
  random_emb <- test[ids_random_emb,]
  all_random_dists <- distance_function(random_emb)
  
  ttest <- t.test(x = all_sampled_dists$value,all_random_dists$value)
  multi2 <- function (x, col = palette(), lwd = 1, lty = 1, xlim, ylim,ylab = "Density", 
                      ...) 
  {
    if (missing(xlim)) {
      xvals <- unlist(lapply(x, function(z) {
        range(z[is.finite(z)], na.rm = TRUE)
      }))
      xlim <- range(xvals[is.finite(xvals)])
    }
    dens <- lapply(x, function(x) {
      density(x[is.finite(x)])
    })
    yvals <- unlist(lapply(dens, function(x) {
      x$y
    }))
    plot(0, type = "n", xlim = xlim, ylim = ylim, ylab = ylab, 
         ...)
    out <- mapply(dens, rep(col, length = length(dens)), rep(lwd, 
                                                             length = length(dens)), rep(lty, length = length(dens)), 
                  FUN = function(x, col, lwd, lty) {
                    lines(x, col = col, lwd = lwd, lty = lty)
                  })
  }
  a <- max(c(density(all_sampled_dists$value)$y,density(all_random_dists$value)$y))
  dir.create(paste0(output_dir,"/task1"),recursive = T)
  png(file=paste0(output_dir,"/task1/","task1_same_vs_random_sig_id.png"),width=7,height=6,units = "in",res=300)
  multi2(list(all_sampled_dists$value,all_random_dists$value),xlab = "Distance", xaxs="i",yaxs="i",ylim = c(0,a+0.01))
  title("A", adj = 0)
  legend("topright", 
         legend = c("graphs from the same sig_id", "graphs from dif sig_id"), 
         col = c('black', 
                 'red'),
         lty = c(1,1),
         bty = "o", 
         pt.cex = 1.5, 
         cex = 0.8, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.01, 0.01))
  dev.off()
  output_task1 <- data.frame(matrix(0,nrow = 1, ncol = 6))
  colnames(output_task1) <- c("mean_same","sd_same","n_same","mean_random","sd_random","n_random")
  output_task1$mean_same <- mean(all_sampled_dists$value,na.rm=T)
  output_task1$sd_same <- sd(all_sampled_dists$value,na.rm=T)
  output_task1$n_same <- nrow(all_sampled_dists)
  output_task1$mean_random <- mean(all_random_dists$value,na.rm=T)
  output_task1$sd_random <- sd(all_random_dists$value,na.rm=T)
  output_task1$n_random <- nrow(all_random_dists)
  write.csv(output_task1,paste0(output_dir,"/task1/","task1_stats.csv"))
  
  # Task 1 TSNE + UMAP
  
  file_info_tsne_1 <- file_info_t1[which(file_info_t1$sig_id %in% sampled_sigs),]
  sampled_sigs <- as.data.frame(sampled_sigs)
  sampled_sigs$label <- c(rep(1,90),2:11)
  file_info_tsne_1 <- left_join(file_info_tsne_1,sampled_sigs, by = c("sig_id"="sampled_sigs"))
  file_info_labels <- file_info_tsne_1 %>% dplyr::select(emb,label)
  emb_tsne <- test[which(test$emb %in% file_info_labels$emb),]
  emb_tsne <- left_join(emb_tsne,file_info_labels,by=c("emb"="emb"))
  emb_tsne$label <- as.factor(emb_tsne$label)
  tsne <- Rtsne(emb_tsne[,-c(1,ncol(emb_tsne))], dims = 2, perplexity = 50, verbose=TRUE, max_iter = 1000,check_duplicates = F)
  df <- data.frame(V1 = tsne$Y[,1], V2 =tsne$Y[,2], label = emb_tsne$label)
  
  png(file=paste0(output_dir,"/task1/","task1_same_sig_id_tsne.png"),width=9,height=9,units = "in",res=300)
  gg1 <- ggplot(df, aes(V1, V2))+
     geom_point(aes(color = label))+aes(group=rev(label))
  print(gg1)
  dev.off()
  
  map <- umap(emb_tsne[,-c(1,ncol(emb_tsne))])
  df_map <- data.frame(V1 = map$layout[,1], V2 = map$layout[,2], label = emb_tsne$label)
  png(file=paste0(output_dir,"/task1/","task1_same_sig_id_umap.png"),width=9,height=9,units = "in",res=300)
  gg_map <- ggplot(df_map, aes(V1, V2))+
    geom_point(aes(color = label))+aes(group=rev(label))
  print(gg_map)
  dev.off()
  print("TASK 1 FINISHED")
  
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
        dup_emb <- slice_sample(filt_embs,1)
        dup_embs <- rbind(dup_emb,dup_embs)
      }
      if (nrow(dup_embs)>1) {
        dup_distance <- distance_function(dup_embs)
        dup_distances <- rbind(dup_distance,dup_distances)
      }
      
    }
  }
  
  sigs <- unique(as.character(file_info$sig_id))
  
  sigs_random <- sample(sigs,35)
  
  random_embs <- NULL
  for (i in 1:length(sigs_random)) {
    filt <- file_info %>% filter(sig_id == sigs_random[i])
    random_emb <- sample_n(test[which(test$emb %in% filt$emb),],1)
    random_embs <- rbind(random_emb,random_embs)
  }
  random_distances <- distance_function(random_embs)
  
  # write results task 2
  b <- max(c(density(dup_distances$value)$y,density(random_distances$value)$y))
  dir.create(paste0(output_dir,"/task2"),recursive = T)
  png(file=paste0(output_dir,"/task2/","task2_duplicate_vs_random_sig_id.png"),width=7,height=6,units = "in",res=300)
  multi2(list(dup_distances$value,random_distances$value),xlab = "Distance", xaxs="i",yaxs="i",ylim = c(0,b+0.01))
  title("B", adj = 0)
  legend("topright", 
         legend = c("graphs from duplicate perts", "random graphs"), 
         col = c('black', 
                 'red'),
         lty = c(1,1),
         bty = "o", 
         pt.cex = 1.5, 
         cex = 0.8, 
         text.col = "black", 
         horiz = F , 
         inset = c(0.01, 0.01))
  dev.off()
  output_task2 <- data.frame(matrix(0,nrow = 1, ncol = 6))
  colnames(output_task2) <- c("mean_same","sd_same","n_same","mean_random","sd_random","n_random")
  output_task2$mean_same <- mean(dup_distances$value,na.rm=T)
  output_task2$sd_same <- sd(dup_distances$value,na.rm=T)
  output_task2$n_same <- nrow(dup_distances)
  output_task2$mean_random <- mean(random_distances$value,na.rm=T)
  output_task2$sd_random <- sd(random_distances$value,na.rm=T)
  output_task2$n_random <- nrow(random_distances)
  write.csv(output_task2,paste0(output_dir,"/task2/","task2_stats.csv"))
  # Task 2 TSNE
  task2_tsne <- file_info_dups[which(file_info_dups$identifier %in% identifier),]
  task2_tsne <- task2_tsne[which(as.character(task2_tsne$emb) %in% as.character(test$emb)),]
  task2_tsne <- task2_tsne %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
  task2_tsne_labels <- task2_tsne %>% dplyr::select(emb,identifier)
  emb_task2_tsne <- test[which(test$emb %in% task2_tsne$emb),]
  emb_task2_tsne <- left_join(emb_task2_tsne,task2_tsne_labels)
  emb_task2_tsne$identifier <- as.factor(emb_task2_tsne$identifier)
  
  tsne2 <- Rtsne(emb_task2_tsne[,-c(1,ncol(emb_task2_tsne))], dims = 2, perplexity=5, verbose=TRUE, max_iter = 1000,check_duplicates = F)
  g <- NULL
  for (i in 1:25) {
    labels <- emb_task2_tsne$identifier
    labels <- as.data.frame(labels)
    labels$label <- "non duplicate"
    labels$label[which((labels$labels %in% identifier[i]))] <- "duplicate" 
    labels$label <- factor(labels$label,levels=c("non duplicate","duplicate"))
    df2 <- data.frame(V1 = tsne2$Y[,1], V2 =tsne2$Y[,2], label = labels$label)
    g2 <- ggplot(df2 %>% mutate(label == "non duplicate") %>% arrange(label), aes(V1, V2))+
      geom_point(aes(color = label),show.legend = T)
    g[[i]] <- g2
  }
  png(file=paste0(output_dir,"/task2/","task_duplicate_sig_id_tnse.png"),width=12,height=9,units = "in",res=300)
  gridExtra::grid.arrange(grobs = g,nrow = 5)
  dev.off()
  
  # UMAP task 2
  umap2 <- umap(emb_task2_tsne[,-c(1,ncol(emb_task2_tsne))])
  gmap <- NULL
  for (i in 1:25) {
    labels <- emb_task2_tsne$identifier
    labels <- as.data.frame(labels)
    labels$label <- "duplicate"
    labels$label[which(!(labels$labels %in% identifier[i]))] <- "non duplicate" 
    labels$label <- factor(labels$label,levels=c("non duplicate","duplicate"))
    df2 <- data.frame(V1 = umap2$layout[,1], V2 =umap2$layout[,2], label = labels$label)
    g2map <- ggplot(df2 %>% mutate(label == "duplicate") %>% arrange(label), aes(V1, V2))+
      geom_point(aes(color = label),show.legend = T)
    gmap[[i]] <- g2map
  }
  png(file=paste0(output_dir,"/task2/","task_duplicate_sig_id_umap.png"),width=12,height=9,units = "in",res=300)
  gridExtra::grid.arrange(grobs = gmap,nrow = 5)
  dev.off()
  #task 2 utility of embeddings compare with gene distance
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
  
  distance_scores <- function(num_table, threshold_count, names) {
    library(GeneExpressionSignature)
    library(tidyverse)
    
    ### rank the table
    table_ranked <- apply(X = -num_table, MARGIN = 2, FUN = rank, ties.method = "random")
    
    ### create the phenodata
    pheno2 <- as.data.frame(colnames(num_table))
    rownames(pheno2) <- colnames(num_table)
    pheno_new <- new("AnnotatedDataFrame",data=pheno2)
    ### create expression set
    expr_set <- new("ExpressionSet",exprs = table_ranked, phenoData=pheno_new)
    ### calculate distances
    distances <- ScoreGSEA(expr_set , threshold_count,"avg")
    colnames(distances) <- names
    rownames(distances) <- names
    return(distances)
  }
  if (genes) {
    ### duplicate gene distances
    thresh <- c(10,15,20,25,30)
    identifier <- unique(as.character(file_info_dups$identifier))
    dup_gene_distances <- readRDS("data/gene_distances_val/dup_gene_distances_for_eval.rds")
    
    # calculate random gene distances
    sigs_random <- as.data.frame(sigs_random)
    sigs_random <- left_join(sigs_random,sig_mapping,by = c("sigs_random"="sig_id2"))
    dist_genes_random <- list(0)
    genes_random <- get_cmap_signatures(sig_ids = sigs_random$sig_id,cmap_path_to_gctx = ds_path,landmark_df = landmark,landmark = TRUE)
    for (j in 1:length(thresh)) {
      dist2 <- distance_scores(num_table = genes_random,threshold_count = thresh[j],names = colnames(genes_random))
      dist_genes_random[[j]] <- dist2
    }
    sam <- 0
    for (k in 1:length(dist_genes_random)) {
      sam <- sam + dist_genes_random[[k]]
    }
    sam <- sam/length(dist_genes_random)
    
    sam[lower.tri(sam,diag = T)] <- 666
    rownames(sam) <- as.character(sigs_random$sigs_random)
    colnames(sam) <- as.character(sigs_random$sigs_random)
    sam <- reshape2::melt(sam)
    sam <- sam %>% filter(value != 666)
    sam$value <- sam$value/2
    
    # write results task 2 genes
    c <- max(c(density(sam$value)$y,density(dup_gene_distances$value)$y))
    png(file=paste0(output_dir,"/task2/","task2_duplicate_vs_random_sig_id_genes.png"),width=7,height=6,units = "in",res=300)
    multi2(list(dup_gene_distances$value,sam$value),xlab = "Distance", xaxs="i",yaxs="i",ylim = c(0,c+0.01))
    title("C", adj = 0)
    legend("topright", 
           legend = c("genes from duplicate perts", "random genes"), 
           col = c('black', 
                   'red'),
           lty = c(1,1),
           bty = "o", 
           pt.cex = 1.5, 
           cex = 0.8, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.01, 0.01))
    dev.off()
    output_task2genes <- data.frame(matrix(0,nrow = 1, ncol = 6))
    colnames(output_task2genes) <- c("mean_same","sd_same","n_same","mean_random","sd_random","n_random")
    output_task2genes$mean_same <- mean(dup_gene_distances$value,na.rm=T)
    output_task2genes$sd_same <- sd(dup_gene_distances$value,na.rm=T)
    output_task2genes$n_same <- nrow(dup_gene_distances)
    output_task2genes$mean_random <- mean(sam$value,na.rm=T)
    output_task2genes$sd_random <- sd(sam$value,na.rm=T)
    output_task2genes$n_random <- nrow(sam)
    write.csv(output_task2genes,paste0(output_dir,"/task2/","task2_stats_genes.csv"))
  }
  
  print("TASK 2 FINISHED")
  # Task 3 similar chemical structure similar embeddings and comparison with gene dist
  if (task3) {
  file_info <- left_join(file_info,sig_mapping,by = c("sig_id"="sig_id2"))
  cell_info <- file_info %>% group_by(cell_id) %>% summarise(count = n_distinct(sig_id))
  cell_info <- cell_info[order(cell_info$count,decreasing = T),]
  if (distance_type == "cosine") {
    # cosine is normalized to 0-1 universally
    distance_function2 <- function(df){
      cos <- cosine(t(df[,-1]),y = NULL)
      cos <- (cos + 1)/(2)
      cos_dist <- 1-cos
      #cos_dist[lower.tri(cos_dist,diag = T)] <- 666
      rownames(cos_dist) <- as.character(df$emb)
      colnames(cos_dist) <- as.character(df$emb)
      cos_dist <- reshape2::melt(cos_dist)
      #cos_dist <- cos_dist %>% filter(value != 666)
      return(cos_dist)
    }
  }
  
  if (distance_type == "euclidian") {
    #this is unnormalized
    distance_function2 <- function(df){
      eu_dist <- dist(df[,-1])
      eu_dist <- as.matrix(eu_dist)
      #eu_dist[lower.tri(eu_dist,diag = T)] <- 666
      rownames(eu_dist) <- as.character(df$emb)
      colnames(eu_dist) <- as.character(df$emb)
      eu_dist <- reshape2::melt(eu_dist)
      #eu_dist <- eu_dist %>% filter(value != 666)
      return(eu_dist)
    }
  }
  dir.create(paste0(output_dir,"/task3"),recursive = T)
  file_info$cell_id <- as.character(file_info$cell_id)
  cell_info$cell_id <-  as.character(cell_info$cell_id)
  for (i in 1:4) {
    file_info_t3 <- file_info %>% filter(cell_id==cell_info$cell_id[i]) %>% filter (!is.na(rdkit))
    rdkit_sig_id_map <- file_info_t3 %>% dplyr::select(rdkit,sig_id,sig_id.y) %>% filter (!is.na(rdkit)) %>% unique()
    rdkit_test <- file_info_t3 %>% dplyr::select(rdkit) %>% filter(!is.na(rdkit)) %>% unique()
    ids_rdkit <- which(rownames(sims) %in% rdkit_test$rdkit)
    sims_filt <- sims[ids_rdkit,ids_rdkit]
    sims_filt <- as.matrix(sims_filt)
    rownames(sims_filt) <- colnames(sims_filt)
    sims_filt <- reshape2::melt(sims_filt)
    
    sims_filt <- left_join(sims_filt,rdkit_sig_id_map,by = c("Var1"="rdkit"))
    sims_filt <- left_join(sims_filt,rdkit_sig_id_map,by = c("Var2"="rdkit"))
    colnames(sims_filt) <- c("rdkit_x","rdkit_y","value_ecfp","sig_id_x1","sig_id_x2","sig_id_y1","sig_id_y2")
    sims_filt$value_ecfp <- 1-sims_filt$value_ecfp
    
    sam <- readRDS(paste0("data/gene_distances_val/",as.character(cell_info$cell_id[i]),"_genes.rds"))
    
    dist_t3 <- left_join(sims_filt,sam,by = c("sig_id_x2"="Var1","sig_id_y2"="Var2"))
    
    # now need to add 5 different embedding distances
    
    for (j in 1:5) {
      file_info_t3_embs <- file_info_t3 %>% group_by(sig_id) %>% mutate(sampled_emb = sample(emb,1)) %>% ungroup() %>%
        dplyr::select(sig_id,sig_id.y,sampled_emb) %>% unique()
      embs <- test[which(test$emb %in% file_info_t3_embs$sampled_emb),]
      dist_embs <- distance_function2(df = embs)
      dist_embs <- left_join(dist_embs,file_info_t3_embs,by = c("Var1"="sampled_emb"))
      dist_embs <- left_join(dist_embs,file_info_t3_embs,by = c("Var2"="sampled_emb"))
      colnames(dist_embs) <- c(paste0("emb_x",j),paste0("emb_y",j),paste0("emb_dist",j),"sig_id_x1","sig_id_x2","sig_id_y1","sig_id_y2")
      dist_t3 <- left_join(dist_t3,dist_embs,by = c("sig_id_x1"="sig_id_x1","sig_id_y1"="sig_id_y1"))
    }
    print(paste0("TASK 3 FINISHED BY ",25*i,"%"))
    dist_t3 <- dist_t3 %>% dplyr::select(rdkit_x,rdkit_y,sig_id_x1,sig_id_y1,value_ecfp,value,emb_x1,emb_y1,emb_dist1,emb_x2,emb_y2,emb_dist2,
                                  emb_x3,emb_y3,emb_dist3,emb_x4,emb_y4,emb_dist4,emb_x5,emb_y5,emb_dist5)
    dist_t3 <- dist_t3 %>% filter(rdkit_x != rdkit_y)
    saveRDS(dist_t3,paste0(output_dir,"/task3/","task3_df_",as.character(cell_info$cell_id[i]),"_all_dists.RDS"))
    
    gene_thresh <- c(0.2,0.4,0.6,0.8,1) * (output_task2genes$mean_same+output_task2genes$sd_same)
    emb_thresh <- c(0.2,0.4,0.6,0.8,1) * (output_task2$mean_same+output_task2$sd_same)
    str_thresh <- c(0.1,0.2,0.3)
    
    gene_prec <- matrix(0,nrow = 4,ncol = 6)
    emb_prec <- matrix(0,nrow = 4, ncol = 6)
    
    for (k in 1:length(str_thresh)) {
      for (l in 1:length(emb_thresh)) {
        n <- length(which(dist_t3$value_ecfp<=str_thresh[k]))
        gene_prec[k,l] <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$value<=gene_thresh[l]))/n
        emb_prec_1 <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$emb_dist1<=emb_thresh[l]))/n
        emb_prec_2 <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$emb_dist2<=emb_thresh[l]))/n
        emb_prec_3 <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$emb_dist3<=emb_thresh[l]))/n
        emb_prec_4 <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$emb_dist4<=emb_thresh[l]))/n
        emb_prec_5 <- length(which(dist_t3$value_ecfp<=str_thresh[k] & dist_t3$emb_dist5<=emb_thresh[l]))/n
        emb_prec[k,l] <- max(emb_prec_1,emb_prec_2,emb_prec_3,emb_prec_4,emb_prec_5)
        
        gene_prec[4,l] <- length(which(dist_t3$value<=gene_thresh[l]))
        n_dist1 <- length(which(dist_t3$emb_dist1<=emb_thresh[l]))
        n_dist2 <- length(which(dist_t3$emb_dist2<=emb_thresh[l]))
        n_dist3 <- length(which(dist_t3$emb_dist3<=emb_thresh[l]))
        n_dist4 <- length(which(dist_t3$emb_dist4<=emb_thresh[l]))
        n_dist5 <- length(which(dist_t3$emb_dist5<=emb_thresh[l]))
        emb_prec[4,l] <- max(n_dist1,n_dist2,n_dist3,n_dist4,n_dist5)
      }
      gene_prec[k,6] <- n
      emb_prec[k,6] <- n
    }
    
    colnames(gene_prec) <- c("20%","40%","60%","80%","100%","no_similar_drugs")
    rownames(gene_prec) <- c("str_thresh_01","str_thresh_02","str_thresh_03","number_lower_than_thresh")
    colnames(emb_prec) <- c("20%","40%","60%","80%","100%","no_similar_drugs")
    rownames(emb_prec) <- c("str_thresh_01","str_thresh_02","str_thresh_03","number_lower_than_thresh")
    write.csv(gene_prec,paste0(output_dir,"/task3/","gene_precision_",as.character(cell_info$cell_id[i]),".csv"),row.names = T) 
    write.csv(emb_prec,paste0(output_dir,"/task3/","emb_precision_",as.character(cell_info$cell_id[i]),".csv"),row.names = T) 
    
    
  
      
  }
  }
  
  }
  

# file info
file_info <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/graph_info_df/file_info_nodups.rds")

# dup file

file_info_dups <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/graph_info_df/file_info_dups.rds")



###
distance_type = "cosine"
###
output_dir <- "/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/validation/validation_transformer_infomax/DGI_DV2_512_seqveq_un_l4/"



ds_path <- "../GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
library(tidyverse)
landmark <- read_tsv(file = "/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/cmap/util_files/cmap_landmark_genes.txt")

sig_mapping <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/graph_info_df/sig_mapping.rds")



### create test embedding df
library(data.table)
test <- fread("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/embeddings/deep_graph_infomax/unsupervised/DGI_JSD_512_seqveq_uniform_un_l4.csv", integer64 = 'double')

#test_files <- as.character(test$X)
#test_files <- as.data.frame(test_files)
#test_files <- left_join(test_files,file_info,by = c("test_files"="files_combined"))
#test_files <- test_files %>% dplyr::select(test_files,emb)
#test_files <- left_join(test_files,file_info_dups,by = c("test_files"="files_combined"))
#test_files <- test_files %>% dplyr::select(test_files,emb.x,emb.y) %>% mutate(emb = if_else(condition = is.na(emb.x),true = emb.y,false = emb.x))
#test_files <- test_files %>% dplyr::select(test_files,emb)
#test <- left_join(test_files,test,by= c("test_files"="X"))

#alevizos nans####
#test <- test %>% filter(!is.na(X0))
#test_names <- read.csv("../../Downloads/val_names.csv")
#test_names <- test_names[,-1]
#test_names <- as.character(test_names)
#test$emb <- test_names
######

test <- test[,-1]
#embs <- read.csv("embeddings/graph2vec/emb_clustered_norm_500.csv")
#test[,1] <- as.character(embs[,1])
colnames(test)[1] <- "emb"
#test <- test %>% mutate(emb = str_remove_all(string = emb,pattern = ".csv"))

# rdkit smiles



sims <- readRDS("/home/rootlocus/Desktop/NTUA/Systems_Biology_Lab/DiplomaThesis/deepSNEM/data/rdkit/rdkit_sims.rds")


eval_emb(test = test, file_info = file_info ,
         distance_type = distance_type,output_dir = output_dir, file_info_dups = file_info_dups, 
         ds_path = NULL,landmark = NULL, sig_mapping = sig_mapping, sims = sims , task3=F, genes = F)
