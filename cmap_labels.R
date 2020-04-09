library(tidyverse)

### cmap drugs

cmap <- readRDS("cmap/pert_id_to_rdkit.rds")

### pert info

pert_info <- read.delim(file = "cmap/GSE92742_Broad_LINCS_pert_info.txt")

pert_info <- pert_info %>%
  filter(pert_type == "trt_cp") %>%
  filter(inchi_key != -666)

### add moa and target from the broad repo

broad_repo <- read.delim(file = "cmap/Repurposing_Hub_export (1).txt" ,skip = 0)

broad_repo <- broad_repo %>%
  mutate(InChIKey = str_split(InChIKey,pattern = ",")) %>% unnest(InChIKey) %>%
  mutate(InChIKey = str_trim(InChIKey)) %>%
  filter(InChIKey != "") %>% unique()

broad_repo <- broad_repo  %>%
  mutate(broad_id = substr(x = as.character(Id),start = 1,stop = 13)) %>%
  mutate(broad_id_old = substr(x = as.character(Deprecated.ID),start = 1,stop = 13)) %>%
  mutate(Name = toupper(Name))

broad_repo <- broad_repo %>%
  mutate(SMILES = str_split(SMILES,pattern = ",")) %>% unnest(SMILES) %>%
  mutate(SMILES = str_trim(SMILES)) %>%
  filter(SMILES != "") %>% unique()

write.csv(broad_repo,"cmap/broad_repo.csv")

broad_repo_rdkit <- read.csv("cmap/broad_repo_rdkit.csv")

broad_repo_rdkit <- broad_repo_rdkit[,-1]
colnames(broad_repo_rdkit) <- c("rdkit","atoms")
broad_repo <- cbind(broad_repo,broad_repo_rdkit)

write.csv(as.character(broad_repo$rdkit),"cmap/repo_rdkit.csv")
write.csv(as.character(cmap$rdkit),"cmap/cmap_smiles.csv")


#### read sims repo to cmap

sims <- read.csv("cmap/cmap_repo_sims.csv")
sims <- sims[,-1]

sims_melt <- reshape::melt(sims)

sims <- t(sims)
max_sims <- apply(sims,1,max,na.rm = T)
length(which(max_sims>=0.6))
