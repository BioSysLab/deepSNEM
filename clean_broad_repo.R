library(tidyverse)

broad_repo <- read.delim2("data/cmap/repo_hub/repurposing_drugs_20200324.txt")
broad_samples <- read.delim2("data/cmap/repo_hub/repurposing_samples_20200324.txt")
broad_repo <- left_join(broad_repo,broad_samples,by = c("ï..pert_iname"="pert_iname"))
broad_repo <- broad_repo %>% select(smiles,moa,target,disease_area,indication) %>% filter (!is.na(smiles)) %>% unique()

# add rdkit to broad repo smiles
repo_rdkit <- read.csv("data/cmap/repo_hub/repo_rdkit_smiles.csv")
repo_rdkit <- repo_rdkit %>% select(X0,Atoms)
colnames(repo_rdkit) <- c("rdkit","atoms")

broad_repo <- cbind(broad_repo,repo_rdkit)
broad_repo <- broad_repo %>% select(rdkit,moa,target,disease_area,indication) %>% filter(!is.na(rdkit)) %>% filter (rdkit!=666) %>% unique()

# remove rdkits without available field
broad_repo <- broad_repo %>% filter(!(moa == "" & target == "" & disease_area == "" & indication == "")) %>% unique()
broad_repo <- broad_repo %>% group_by(rdkit) %>% mutate(count = n()) %>% ungroup()

broad_repo_doubles <- broad_repo %>% filter (count>1)
broad_repo_singles <- broad_repo %>% filter(count == 1)

broad_repo_doubles$keep <- NULL
for (i in 1:nrow(broad_repo_doubles)) {
  
  sum <- sum(broad_repo_doubles$moa[i] == "",broad_repo_doubles$target[i] == "",broad_repo_doubles$disease_area[i] == "",broad_repo_doubles$indication[i] == "")
  broad_repo_doubles$keep[i] <- sum
  
}
broad_repo_doubles <- broad_repo_doubles %>% group_by(rdkit) %>% filter(keep == min(keep)) %>% ungroup()
broad_repo_doubles <- broad_repo_doubles %>% mutate(gg = paste0(moa,target,disease_area,indication)) %>% mutate(gglen = nchar(gg))
broad_repo_doubles <- broad_repo_doubles %>% group_by(rdkit) %>% filter(gglen == max(gglen))%>% ungroup()
broad_repo_doubles <- broad_repo_doubles %>% select(rdkit,moa,target,disease_area,indication,count)

broad_repo <- bind_rows(broad_repo_singles,broad_repo_doubles) %>% select(rdkit,moa,target,disease_area,indication)

saveRDS(broad_repo,"data/cmap/repo_hub/broad_repurposing.rds")

write.csv(as.character(broad_repo$rdkit),"data/cmap/repo_hub/repo_rdkit.csv")
