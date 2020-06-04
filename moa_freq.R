library(tidyverse)

allpairs <- readRDS("data/graph_additional/pairs/splits/split3/allpairs3.rds")

labeled <- allpairs %>% filter(label==1) %>% 
  filter(!is.na(moa_v1.x)) %>% filter(!is.na(moa_v1.y)) %>% 
  filter(moa_v1.x != "") %>% filter(moa_v1.y != "")

freq_labels <- labeled %>% group_by(moa_v1.x) %>% summarise(count = n()/nrow(labeled)) %>% arrange(desc(count))

saveRDS(freq_labels,"data/graph_additional/pairs/splits/split3/moa_df.rds")
write.csv(freq_labels,"data/graph_additional/pairs/splits/split3/moa_df.csv")
