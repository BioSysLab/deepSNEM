library(tidyverse)

files <- list.files(path = "graphs",recursive = T,full.names = T)
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

files_nets <- files_nets %>% select(files) %>% 
  mutate(weighted = grepl(pattern = paste(match_weight,collapse = "|"),x = files))

files_nets_weighted <- files_nets %>% filter(weighted == T)


files_nets_weighted <-  files_nets_weighted %>% mutate(type = if_else(condition = grepl(pattern = "Attributes",x = files),true = "nodes",false = "edges"))
files_nets_weighted <- files_nets_weighted %>% select(-weighted)

files_nets <- files_nets %>% filter(weighted == F)

files_nets <- files_nets %>% select(files) %>%
  mutate(type = if_else(condition = grepl(pattern = "nodesActivity",x = files),true = "nodes",false = "edges"))

# from both extract sig id


files_nets_weighted <- files_nets_weighted %>%
  mutate(sig_id = qdapRegex::ex_between(files, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 2)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))

files_nets <- files_nets %>%
  mutate(sig_id = qdapRegex::ex_between(files, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 2)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))


no_weights <- files_nets %>%  select(files,sig_id) %>% group_by(sig_id) %>%
  summarise(count = n()) %>% arrange(count) %>% mutate(count = count/2)

files_nets <- left_join(files_nets,no_weights)

weights <- files_nets_weighted %>% select(files,sig_id) %>% group_by(sig_id) %>%
  summarise(count = n()) %>% arrange(count) %>% mutate(count = count/2)

files_nets_weighted <- left_join(files_nets_weighted,weights)

### create one file

sig_ids <- files_nets %>% select(sig_id,count) %>% unique() 

for (i in 1:nrow(sig_ids)) {
  
  filt <- files_nets %>% filter(sig_id == sig_ids$sig_id[i])
  filt_nodes <- filt %>% filter(type == "nodes") %>% select(-type) 
  colnames(filt_nodes) <- c("files_nodes","sig_id","count")
  filt_edges <- filt %>% filter(type == "edges") %>% select(-type)
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
      
      dir.create(paste0("graphs_combined/",as.character(sig_ids$sig_id[i])))
      write.csv(combo,paste0("graphs_combined/",as.character(sig_ids$sig_id[i]),"/graph_",as.character(filt_nodes$model_no[j]),".csv"))
      
    }
  }
  if (all(filt_nodes$model_no != filt_edges$model_no)) {
    print("nope")
  }
  print(i)
}

files_combined <- list.files(path = "graphs_combined",recursive = T,full.names = T)
files_combined <- as.data.frame(files_combined)
files_combined$files_combined <- as.character(files_combined$files_combined)
files_combined <- files_combined %>%
  mutate(sig_id = qdapRegex::ex_between(files_combined, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 1)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))

number <- files_combined %>% group_by(sig_id) %>% summarize(count = n()) %>% arrange()

files_combined <- left_join(files_combined,number)

sigs <- read.csv("SNAC_database.csv")
sigs <- sigs[,-1]
sigs <- sigs %>% mutate(sig_id = str_replace_all(string = sig_id,pattern = ":",replacement = "_"))
files_combined <- left_join(files_combined,sigs, by = "sig_id")


file_info <- readRDS("data/file_info.rds")

files_info <- files_info %>% mutate(emb = files_combined) %>%
  mutate(emb = str_replace_all(string = emb,pattern = "graphs_combined/",replacement = "")) 

files_info <- files_info %>%
  mutate(emb = str_replace_all(string = emb, pattern = "/graph",replacement = "_emb"))

files_info <- files_info %>%
  mutate(emb = str_replace_all(string = emb, pattern = ".csv",replacement = "")) %>% mutate(emb = str_trim(emb))

saveRDS(files_info,"file_info.rds")
write.csv(files_info,"file_info.csv")
