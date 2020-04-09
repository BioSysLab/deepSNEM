library(tidyverse)
sig <- read.delim("C:/Users/user/Documents/MYC/data/cmap_data/GSE92742_Broad_LINCS_sig_info.txt")
sig_metrics <- read.delim("C:/Users/user/Documents/MYC/data/cmap_data/GSE92742_Broad_LINCS_sig_metrics.txt")  

cell <- sig  %>%
  filter(pert_type == "trt_cp") %>%
  group_by(pert_iname,cell_id) %>%
  mutate(count = n_distinct(sig_id)) %>%
  ungroup() %>% filter(count>1)

cell <- left_join(cell,sig_metrics)
cell$quality <- 100

cell <- cell %>%
  mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
         quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
         quality = factor(quality))

cell <- cell %>% filter(quality==1)


cell <- cell %>%
  group_by(pert_iname,cell_id) %>%
  mutate(count = n_distinct(sig_id),
         count_dose = n_distinct(pert_dose),
         count_time = n_distinct(pert_time)) %>% 
  ungroup() %>% filter(count > 1) %>%
  filter(count_dose == 1) %>% filter(count_time == 1)

cell <- cell %>%
  group_by(pert_iname,cell_id) %>%
  mutate(count = n_distinct(sig_id)) %>% ungroup() %>%
  mutate(sig_id2 = str_replace_all(string = sig_id,pattern = ":",replacement = "_"))
cell <- cell %>% mutate(identifier = paste0(pert_iname,cell_id,pert_dose,pert_time))

saveRDS(cell,"duplicate_sigs_all.RDS")
cell_torun <- cell[-which(cell$sig_id2 %in% file_info$sig_id),]
saveRDS(cell_torun,"duplicates_for_carnival.RDS")
landmark <- read_tsv(file = "C:/Users/user/Documents/MYC/myc_cmap_pathways/cmap_landmark_genes.txt")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
genes <- get_cmap_signatures(cmap_path_to_gctx = ds_path,sig_ids = as.character(cell_torun$sig_id),landmark = T,landmark_df = landmark)

tfs <- tf_enrichment(sigs = as.character(cell_torun$sig_id),cmap_path = ds_path,landmark = landmark)

saveRDS(tfs,"tfs_dups_carnival.RDS")


cell <- sig  %>%
  filter(pert_type == "trt_cp") 

cell <- left_join(cell,sig_metrics)
cell$quality <- 100

cell <- cell %>%
  mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
         quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
         quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
         quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
         quality = factor(quality))

cell <- cell %>% filter(quality==1)

cell <- cell %>% dplyr::select(sig_id)
cell <- cell %>% mutate(sig_id2 = str_replace_all(string = sig_id,pattern = ":",replacement = "_"))
saveRDS(cell,"sig_mapping.RDS")
