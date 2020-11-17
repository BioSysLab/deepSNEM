library(tidyverse)

df <- readRDS("data/graph_additional/df1.rds")
df2 <- readRDS("data/graph_additional/df2.rds")
df <- bind_rows(df,df2)

df$same_label[df$moa.x==df$moa.y] <- 1

# extract test set

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
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
emb <- read.csv("embeddings/ged_distance/ged_512_cosine.csv")
emb <- emb[,-1]
colnames(emb)[1] <- "emb"
labels <- readRDS("data/cmap/labels/labels_first_pass.rds")
emb_proc <- prepape_embs(emb = emb,file_info = file_info,labels = labels,keep_one = F ,ave = T,n_emb = (ncol(emb)-1))

emb_proc <- emb_proc %>% filter(!is.na(moa_v1))

moa <- emb_proc %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_vis <- moa$moa_v1[1:10]
emb_proc <- emb_proc[which(as.character(emb_proc$moa_v1) %in% moa_vis),]

emb_drugs <- emb_proc %>% group_by(rdkit) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
emb_drugs$id <- 1:nrow(emb_drugs)
emb_drugs <- left_join(emb_drugs,labels,by=c("rdkit"="rdkit_graph"))
emb_id <- emb_drugs %>% select(rdkit,id) %>% unique()
emb_size <- 512
perpl_emb <- 3
iter <- 2000
init_dim <- 80
tsne_emb <- Rtsne(scale(emb_proc[,2:(emb_size+1)]), dims = 2, perplexity=perpl_emb, verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
emb_proc <- left_join(emb_proc,emb_id)
labels <- as.character(emb_proc$id)
#labels[which(emb_proc$rdkit == "COC(=O)Nc1ccc(-c2nc(N3CCOCC3)c3cnn(C4CCN(C(=O)OC)CC4)c3n2)cc1")] <- "x"
df_emb <- data.frame(V1 = tsne_emb$Y[,1], V2 =tsne_emb$Y[,2], label = as.factor(emb_proc$moa_v1), name = labels)
gtsne <- ggplot(df_emb, aes(V1, V2))+
  geom_point(aes(color = label),show.legend = T) + scale_color_discrete()+
  geom_text(aes(label=name),hjust=0, vjust=0,size = 0.5)
png(file=paste0("vis/test_all.png"),width=9,height=9,units = "in",res=600)
print(addSmallLegend(gtsne))
dev.off()


test_rdkit <- c("O=C(CCCCCCC(=O)Nc1ccccc1)NO","CCc1c2c(nc3ccc(O)cc13)-c1cc3c(c(=O)n1C2)COC(=O)C3(O)CC",
                "C[C@@H]1O[C@@H](O[C@H]2C[C@@H](O)[C@]3(CO)C4C(CC[C@]3(O)C2)[C@@]2(O)CC[C@H](C3=CC(=O)OC3)[C@@]2(C)C[C@H]4O)[C@H](O)[C@H](O)[C@H]1O",
                "COc1ccc(-c2ccc3c(N4CCOC[C@@H]4C)nc(N4CCOC[C@@H]4C)nc3n2)cc1CO",
                "CC[C@H]1CN2CCc3cc(OC)c(OC)cc3C2C[C@@H]1CC1NCCc2cc(OC)c(OC)cc21",
                "CC1CCN(CCCC(=O)c2ccc(F)cc2)CC1",as.character(emb_id$rdkit[emb_id$id==10]))

test_rdkit <- as.data.frame(test_rdkit)
rdkit_sig <- emb_proc %>% select(sig_id,rdkit) %>% unique()

test_rdkit <- left_join(test_rdkit,rdkit_sig,by = c("test_rdkit"="rdkit"))

id_test <- unique(c(which(df$sig_id.x %in% test_rdkit$sig_id),which(df$sig_id.y %in% test_rdkit$sig_id)))

df_train <- df[-id_test,]
df_test <- df[id_test,]

sig_id_rdkit <- file_info %>% select(sig_id,rdkit) %>% unique()

df_train <- left_join(df_train,sig_id_rdkit,by=c("sig_id.x"="sig_id"))
df_train <- left_join(df_train,sig_id_rdkit,by=c("sig_id.y"="sig_id"))
rdkit_moa_v1 <- labels %>% select(rdkit_graph,moa_v1)

df_train <- left_join(df_train,rdkit_moa_v1,by=c("rdkit.x"="rdkit_graph"))
df_train <- left_join(df_train,rdkit_moa_v1,by=c("rdkit.y"="rdkit_graph"))
df_train$label <- 0
df_train$label[df_train$moa_v1.x==df_train$moa_v1.y] <- 1
df_train$label[df_train$identifier.x==df_train$identifier.y] <- 1
same_label <- df_train %>% filter(label==1)

saveRDS(df_train,"data/graph_additional/df_train.rds")
saveRDS(df_test,"data/graph_additional/df_test.rds")
saveRDS(test_rdkit,"data/graph_additional/test_sig.rds")
write.csv(df_train,"data/graph_additional/df_train.csv")
write.csv(df_test,"data/graph_additional/df_test.csv")
write.csv(test_rdkit,"data/graph_additional/test_sig.csv")


labels <- readRDS("data/cmap/labels/labels_first_pass.rds")


length(unique(c(df_train$rdkit.x,df_train$rdkit.y)))
length(unique(c(df_train$sig_id.x,df_train$sig_id.y)))
length(unique(c(df_train$graph.x,df_train$graph.y)))

length(unique(c(df$sig_id.x,df$sig_id.y)))
