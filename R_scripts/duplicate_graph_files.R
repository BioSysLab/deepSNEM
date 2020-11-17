library(tidyverse)

dups <- readRDS("duplicate_sigs_all.RDS")




#file_info <- readRDS("data/file_info.rds")

files <- list.files(path = "dublicates",recursive = T,full.names = T)
files <- as.data.frame(files)

files <- files %>%
  mutate(meas = grepl(pattern = "meas_",x = files),
         log = grepl(pattern = ".log",x = files),
         time = grepl(pattern = "elapsed_time.txt",x = files),
         res = grepl(pattern = "results_CARNIVAL.Rdata",x = files),
         empty = grepl(pattern = "emptyNetwork",x = files))

files_nets <-  files %>% filter(meas == F & log == F & time == F & res == F & empty == F)

#
match_weight <- c("nodesAttributes_1.txt","weightedModel_1.txt")

files_nets <- files_nets %>% dplyr::select(files) %>% 
  mutate(weighted = grepl(pattern = paste(match_weight,collapse = "|"),x = files))

files_nets <- files_nets %>% filter(weighted == F)

files_nets <- files_nets %>% dplyr::select(files) %>%
  mutate(type = if_else(condition = grepl(pattern = "nodesActivity",x = files),true = "nodes",false = "edges"))

files_nets <- files_nets %>%
  mutate(sig_id = qdapRegex::ex_between(files, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 2)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))


no_weights <- files_nets %>%  dplyr::select(files,sig_id) %>% group_by(sig_id) %>%
  summarise(count = n()) %>% arrange(count) %>% mutate(count = count/2)

files_nets <- left_join(files_nets,no_weights)
files_nets <- files_nets %>% filter(!is.na(sig_id))
### create one file

sig_ids <- files_nets %>% dplyr::select(sig_id,count) %>% unique() 
library(dplyr)
for (i in 1:nrow(sig_ids)) {
  
  filt <- files_nets %>% filter(sig_id == sig_ids$sig_id[i])
  filt_nodes <- filt %>% filter(type == "nodes") %>% dplyr::select(-type) 
  colnames(filt_nodes) <- c("files_nodes","sig_id","count")
  filt_edges <- filt %>% filter(type == "edges") %>% dplyr::select(-type)
  colnames(filt_edges) <- c("files_edges","sig_id","count")
  
  filt_nodes <- filt_nodes %>% mutate(model_no = qdapRegex::ex_between(files_nodes, "model", ".txt")) %>%
    mutate(model_no = sapply(model_no, "[", 1)) %>% unique() 
  filt_edges <- filt_edges %>% mutate(model_no = qdapRegex::ex_between(files_edges, "model", ".tsv")) %>%
    mutate(model_no = sapply(model_no, "[", 1)) %>% unique() 
  if (all(filt_nodes$model_no == filt_edges$model_no)) {
    for (j in 1:nrow(filt_nodes)) {
      
      nodes <- read.delim(as.character(filt_nodes$files_nodes[j]))
      nodes <- nodes %>% mutate(Activity = round(Activity))
      
      edges <- read.delim(as.character(filt_edges$files_edges[j]))
      
      combo <- left_join(edges,nodes,by = c("Node1"="Nodes"))
      
      combo <- left_join(combo,nodes,by = c("Node2"="Nodes"))
      
      colnames(combo) <- c("node1","sign","node2","activity1","activity2")
      
      dir.create(paste0("dupl/",as.character(sig_ids$sig_id[i])))
      write.csv(combo,paste0("dupl/",as.character(sig_ids$sig_id[i]),"/graph_",as.character(filt_nodes$model_no[j]),".csv"))
      
    }
  }
  if (all(filt_nodes$model_no != filt_edges$model_no)) {
    print("nope")
  }
  print(i)
}


files_combined <- list.files(path = "dupl",recursive = T,full.names = T)
files_combined <- as.data.frame(files_combined)
files_combined$files_combined <- as.character(files_combined$files_combined)
files_combined <- files_combined %>%
  mutate(sig_id = qdapRegex::ex_between(files_combined, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 1)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))

number <- files_combined %>% group_by(sig_id) %>% summarize(count = n()) %>% arrange()

files_combined <- left_join(files_combined,number)
files_combined <- files_combined %>% mutate(files_combined = str_replace_all(string = files_combined,pattern = "dupl/",replacement = "graphs_combined/"))

sigs <- read.csv("data/SNAC_database.csv")
sigs <- sigs[,-1]
sigs <- sigs %>% mutate(sig_id = str_replace_all(string = sig_id,pattern = ":",replacement = "_"))
files_combined <- left_join(files_combined,sigs, by = "sig_id")


files_info <- files_combined

files_info <- files_info %>% mutate(emb = files_combined) %>%
  mutate(emb = str_replace_all(string = emb,pattern = "graphs_combined/",replacement = "")) 

files_info <- files_info %>%
  mutate(emb = str_replace_all(string = emb, pattern = "/graph",replacement = "_emb"))

files_info <- files_info %>%
  mutate(emb = str_replace_all(string = emb, pattern = ".csv",replacement = "")) %>% mutate(emb = str_trim(emb))
#saveRDS(files_info,"sig_id_dupl_g2v.rds")
#saveRDS(files_info,"file_info.rds")

file_info_dups <- files_info[which(files_info$sig_id %in% dups$sig_id2),]
file_info_nodups <- files_info[-which(files_info$sig_id %in% dups$sig_id2),]

saveRDS(file_info_nodups,"file_info_nodups.rds")
saveRDS(file_info_dups,"file_info_dups.rds")
