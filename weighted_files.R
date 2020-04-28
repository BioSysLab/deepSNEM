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

files_nets_weighted <- files_nets_weighted %>%
  mutate(sig_id = qdapRegex::ex_between(files, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 2)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))

weights <- files_nets_weighted %>% select(files,sig_id) %>% group_by(sig_id) %>%
  summarise(count = n()) %>% arrange(count) %>% mutate(count = count/2)

files_nets_weighted <- left_join(files_nets_weighted,weights)

# save to one csv file 

sig_ids <- files_nets_weighted %>% select(sig_id,count) %>% unique() 

for (i in 1:nrow(sig_ids)) {
  #i <- 109
  filt <- files_nets_weighted %>% filter(sig_id == sig_ids$sig_id[i])
  filt_nodes <- filt %>% filter(type == "nodes") %>% select(-type) 
  colnames(filt_nodes) <- c("files_nodes","sig_id","count")
  filt_edges <- filt %>% filter(type == "edges") %>% select(-type)
  colnames(filt_edges) <- c("files_edges","sig_id","count")
  
  nodes <- read.delim(as.character(filt_nodes$files_nodes[1]))
  nodes <- nodes %>% select(Node,UpAct,DownAct) %>% filter(!(UpAct == 0 & DownAct == 0 ))
      
  edges <- read.delim(as.character(filt_edges$files_edges[1]))
  edgenodes <- unique(c(as.character(edges$Node1),as.character(edges$Node2)))
  nodes$Node[which(!(nodes$Node %in% edgenodes))]
  if (length(edgenodes)!=nrow(nodes)) {
    print(nodes$Node[which(!(nodes$Node %in% edgenodes))])
    print(i)
  }
  
  combo <- left_join(edges,nodes,by = c("Node1"="Node"))
      
  combo <- left_join(combo,nodes,by = c("Node2"="Node"))
      
  colnames(combo) <- c("node1","sign","node2","weight","upact1","downact1","upact2","downact2")
  combo <- combo %>% mutate(upact1 = upact1/100,
                            downact1 = downact1/100,
                            upact2 = upact2/100,
                            downact2 = downact2/100,
                            weight = weight/100)    
  dir.create(paste0("graphs_weighted/",as.character(sig_ids$sig_id[i])))
  write.csv(combo,paste0("graphs_weighted/",as.character(sig_ids$sig_id[i]),"/graph.csv"))
  #print(i)
}

files_weighted <- list.files(path = "graphs_weighted",recursive = T,full.names = T)
files_weighted <- as.data.frame(files_weighted)
files_weighted$files_weighted <- as.character(files_weighted$files_weighted)
files_weighted <- files_weighted %>%
  mutate(sig_id = qdapRegex::ex_between(files_weighted, "/", "/") ) %>% 
  mutate(sig_id = sapply(sig_id, "[", 1)) %>% unique() %>% mutate(sig_id = str_trim(sig_id))

number <- files_weighted %>% group_by(sig_id) %>% summarize(count = n()) %>% arrange()

files_weighted <- left_join(files_weighted,number)

# read sigs 
sigs <- read.csv("data/SNAC_database.csv")
sigs <- sigs[,-1]
dups <- readRDS("duplicate_sigs_all.RDS")

pert_to <- readRDS("cmap/util_files/pert_id_to_rdkit.rds")

dups <- left_join(dups,pert_to,by = "pert_id")

sigs <- sigs %>% select(sig_id,rdkit,cell_id,pert_iname,pert_dose,pert_time,pert_id,atoms)
dups <- dups %>% select(sig_id,rdkit,cell_id,pert_iname,pert_dose,pert_time,pert_id,atoms)

all <- rbind(sigs,dups) %>% unique()
all <- all %>% mutate(sig_id = str_replace_all(string = sig_id,pattern = ":",replacement = "_"))
files_weighted <- left_join(files_weighted,all,by = "sig_id")

files_weighted$pert_dose <- as.numeric(files_weighted$pert_dose)
files_weighted <- files_weighted %>% mutate(pert_dose=round(pert_dose,digits = 1)) %>% unique()

files_weighted <- files_weighted %>% mutate(emb = files_weighted) %>%
  mutate(emb = str_replace_all(string = emb,pattern = "graphs_weighted/",replacement = "")) 

files_weighted <- files_weighted %>%
  mutate(emb = str_replace_all(string = emb, pattern = "/graph",replacement = "_emb"))

files_weighted <- files_weighted %>%
  mutate(emb = str_replace_all(string = emb, pattern = ".csv",replacement = "")) %>% mutate(emb = str_trim(emb))

files_weighted <- files_weighted %>% mutate(identifier = paste0(pert_iname,cell_id,pert_dose,pert_time))
sig_map <- readRDS("sig_mapping.RDS")
colnames(sig_map) <- c("sig_id_original","sig_id")
files_weighted <- left_join(files_weighted,sig_map,by = "sig_id")
saveRDS(files_weighted,"file_info_weighted.rds")
write.csv(files_weighted,"file_info_weighted.csv")
