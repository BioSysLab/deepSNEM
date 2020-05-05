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
