length(unique(emb_proc$sig_id))

emb_proc <- emb_proc %>% filter(!is.na(moa_v1)) %>% filter(moa_v1 != "")

emb_proc <- emb_proc %>% group_by(moa_v1) %>% mutate(count_moa = n_distinct(sig_id)) %>% ungroup()
emb_proc <- emb_proc %>% filter(count_moa >= 3)

emb_proc_df <- emb_proc %>% select(sig_id_original,rdkit,rdkit_broad,moa,target,disease_area,indication,moa_v1,count) 

emb_proc <- emb_proc[,1:176]

write.csv(emb_proc_df,"../deepSNEM_personal/emb_proc_df_175_tfs_all.csv",row.names = F)
write.csv(emb_proc,"../deepSNEM_personal/emb_proc_175_tfs_all.csv",row.names = F)

file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
labels <- readRDS("data/cmap/labels/labels_first_pass.rds")
allpairs <- readRDS("embeddings/ged_distance_semi/split3/allpairs3.rds")

samples <- read.csv("../deepsnem_drive/data/samples_all.csv")
samples <- samples[,-1]
samples <- as.data.frame(samples)
samples <- left_join(samples,file_info, by =c("samples"="files_combined"))

samples_dups <- samples %>% filter(is.na(sig_id))
samples_nodups <- anti_join(samples,samples_dups)

samples_dups <- samples_dups %>% dplyr::select(samples)
samples_dups <- left_join(samples_dups,file_info_dups,by=c("samples"="files_combined"))

rdkit_map <- readRDS("data/cmap/util_files/pert_id_to_rdkit.rds")
samples_dups <- left_join(samples_dups,rdkit_map,by=c("pert_id"="pert_id"))
samples_dups <- samples_dups %>% filter(!is.na(rdkit)) %>% filter(rdkit != "")

samples_nodups <- samples_nodups %>% select(samples,sig_id,emb,rdkit)
samples_dups <- samples_dups %>% select(samples,sig_id,emb,rdkit)

samples <- bind_rows(samples_dups,samples_nodups)

samples <- left_join(samples,labels,by =c("rdkit"="rdkit_graph"))

samples <- samples %>% group_by(moa_v1) %>% mutate(count = n_distinct(sig_id)) %>% ungroup()

write.csv(samples,"data/samples_with_labels_g2v_128.csv",row.names = F)


df <- read.csv("data/ged_unsupervised/sig_emb_df_ged.csv")

df <- df %>% select(-count)

length(unique(df$sig_id))
length(which(df$sig_id %in% file_info$sig_id))
info <- file_info %>% select(sig_id,files_combined,emb) %>% unique()

df <- left_join(df,info,by=c("sig_id"="sig_id"))

length(unique(df$sig_id))
length(which(df$moa_v1==""))
length(which(is.na(df$moa_v1)))

write.csv(df,"data/graph_classification/graph_classification_all.csv",row.names = F)
