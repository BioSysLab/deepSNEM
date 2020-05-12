library(tidyverse)

### cmap drugs

cmap <- readRDS("data/cmap/util_files/pert_id_to_rdkit.rds")


# read broad repurposing data

broad_repo <- readRDS("data/cmap/repo_hub/broad_repurposing.rds")



#### read the quality 1 smiles for which we hve the graphs

data_dups <- readRDS("data/graph_info_df/file_info_dups.rds")
data <- readRDS("data/graph_info_df/file_info_nodups.rds")

data <- data %>% select(rdkit) %>% filter(!is.na(rdkit)) %>% unique()
data_dups <- left_join(data_dups,cmap, by = "pert_id")
data_dups <- data_dups %>% select(rdkit) %>% filter(!is.na(rdkit)) %>% unique()

### rdkit with graph available
q1_rdkits <- bind_rows(data,data_dups) %>% unique()

sims <- read.csv("data/cmap/repo_hub/deepsnem_graph_repo_sims.csv") 
sims <- sims[,-1]
sims <- sims > 0.99
sims <- sims+0
rownames(sims) <- as.character(q1_rdkits$rdkit)
sims <- sims[which(rowSums(sims)!=0),]

# add labels to rdkit wiith graphs
labels <- data.frame(matrix(0,nrow=nrow(sims),ncol=ncol(broad_repo)+1))
colnames(labels) <- c("rdkit_graph","rdkit_broad","moa","target","disease_area","indication")
broad_repo$rdkit <- as.character(broad_repo$rdkit)
broad_repo$moa <- as.character(broad_repo$moa)
broad_repo$target <- as.character(broad_repo$target)
broad_repo$disease_area <- as.character(broad_repo$disease_area)
broad_repo$indication <- as.character(broad_repo$indication)
for (i in 1:nrow(sims)) {
  id <- which(sims[i,] == 1)[1]
  labels[i,"rdkit_graph"] <- rownames(sims)[i]
  labels[i,2:ncol(labels)] <- as.character(broad_repo[id,])
}

saveRDS(labels,"data/cmap/labels.rds")
write.csv(labels,"data/cmap/labels.csv")
labels <- readRDS("data/cmap/labels.rds")
moa <- labels %>% group_by(moa) %>% summarise(count = n()) %>% arrange(desc(count))
target <- labels %>% group_by(target) %>% summarise(count = n()) %>% arrange(desc(count))
disease_area <- labels %>% group_by(disease_area) %>% summarise(count = n()) %>% arrange(desc(count))
indication <- labels %>% group_by(indication) %>% summarise(count = n()) %>% arrange(desc(count))

# work on moa grouping

dna_rna_damage <- c("topoisomerase inhibitor","RNA synthesis inhibitor|topoisomerase inhibitor",
                              "chelating agent|topoisomerase inhibitor","RNA synthesis inhibitor","DNA synthesis inhibitor",
                              "DNA alkylating agent","DNA alkylating agent|DNA synthesis inhibitor","DNA inhibitor",
                              "DNA replication inhibitor|STAT inhibitor")

labels$moa_v1 <- labels$moa
labels$moa_v1[which(labels$moa %in% dna_rna_damage)] <- "dna_rna_damage"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))


# if group together agonist antagonist
# file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")
file_info <- file_info %>% dplyr::select(files_combined,sig_id,rdkit,cell_id,count.x,emb)
# embs
emb <- read.csv("embeddings/graph2vec/emb_activity_1_epoch.csv")
colnames(emb)[1] <- "emb"
#emb <- emb %>% mutate(emb = str_remove_all(string = emb,pattern = ".csv"))
# keep only 1 emb per sig id for clustering
# keep the respective embeddings

emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
file_info <- file_info[which(as.character(file_info$emb) %in% as.character(emb$emb)),]

file_info <- file_info %>% group_by(sig_id) %>% sample_n(1) %>% ungroup()
emb <- emb[which(as.character(emb$emb) %in% as.character(file_info$emb)),]
dopamine <- c("dopamine receptor antagonist","dopamine receptor antagonist|serotonin receptor antagonist")
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
landmark <- read_tsv(file = "data/cmap/util_files/cmap_landmark_genes.txt")
sig_map <- readRDS("data/graph_info_df/sig_mapping.rds")
group_name <- dopamine
output_dir <- paste0("group_check/","dopamine")
group_check(group = dopamine,group_name = "dopamine",original_big_group = "dopamine receptor antagonist",
            labels = labels, emb_test = emb,file_info = file_info,emb_size = 128,
            output_dir = output_dir, label_space_only = F, tsne_perpl = 80, init_dim = 50,umap_n = 50, cell_specific = F,cell_line = "A375",
            genes = T,
            ds_path = ds_path, landmark = landmark, sig_map = sig_map, tsne_perpl_genes = 50, init_dim_genes = 50)

labels$moa_v1[which(labels$moa_v1 %in% dopamine)] <- "dopamine_antagonist"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))


#cycloxygenase group

cox <- c("cyclooxygenase inhibitor","cyclooxygenase inhibitor|FAAH inhibitor|TRPV antagonist",
         "cyclooxygenase inhibitor|lipoxygenase inhibitor|prostanoid receptor antagonist",
         "cyclooxygenase inhibitor|PPAR receptor agonist","cyclooxygenase inhibitor|prostanoid receptor agonist")

group_name <- "cox"
output_dir <- paste0("group_check/","cox")
labels <- readRDS("data/cmap/labels.rds")
group_check(group = cox,group_name = group_name,original_big_group = "cyclooxygenase inhibitor",
            labels = labels, emb_test = emb,file_info = file_info,emb_size = 128,
            output_dir = output_dir, label_space_only = F, tsne_perpl = 80, init_dim = 50,umap_n = 50, cell_specific = F,cell_line = "A375",
            genes = T,
            ds_path = ds_path, landmark = landmark, sig_map = sig_map, tsne_perpl_genes = 50, init_dim_genes = 50)


labels$moa_v1[which(labels$moa_v1 %in% cox)] <- "cox inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# adrenergic agonist

adrenergic_agonist <- c("adrenergic receptor agonist","adrenergic receptor agonist|imidazoline receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% adrenergic_agonist)] <- "adrenergic_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# ATP synthesis inhibitor

atp_synthesis <- c("ATPase inhibitor","ATP synthase inhibitor|ATPase inhibitor","ATP synthase inhibitor",
                   "ATPase inhibitor|TRPV agonist")
labels$moa_v1[which(labels$moa_v1 %in% atp_synthesis)] <- "atp_synthesis_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# EGFR inhibitor

egfr <- c("EGFR inhibitor","EGFR inhibitor|tyrosine kinase inhibitor","EGFR inhibitor|protein tyrosine kinase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% egfr)] <- "EGFR_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# adrenergic receptor antagonist

adrenergic_antagonist <- c("adrenergic receptor antagonist","adrenergic receptor antagonist|prolactin inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% adrenergic_antagonist)] <- "adrenergic_receptor_antagonist"

# create multiples

id_1 <- which(as.character(labels$moa) == "adrenergic receptor antagonist|serotonin receptor antagonist|serotonin reuptake inhibitor")

multi <- labels[id_1,]
multi <- bind_rows(multi,multi)

multi$moa_v1[1] <- "adrenergic_receptor_antagonist"
multi$moa_v1[2] <- "serotonin receptor antagonist"

labels <- labels[-id_1,]

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# glutamate receptor antagonist

glutamate_antagonist <- c("glutamate receptor antagonist","glutamate inhibitor","glutamate receptor antagonist|glutaminase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% glutamate_antagonist)] <- "glutamate_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# mtor + pi3k + plk

multi2 <- labels[which(labels$moa == "DNA dependent protein kinase inhibitor|mTOR inhibitor|phosphodiesterase inhibitor|PI3K inhibitor|PLK inhibitor"),]

multi2 <- bind_rows(multi2,multi2,multi2)
multi2$moa_v1 <- c("mtor_inhibitor","plk_inhibitor","pi3k_inhibitor")
labels <- labels[-which(labels$moa == "DNA dependent protein kinase inhibitor|mTOR inhibitor|phosphodiesterase inhibitor|PI3K inhibitor|PLK inhibitor"),]
# mtor + pi3k 

multi3 <- labels[which(labels$moa == "DNA dependent protein kinase inhibitor|mTOR inhibitor|PI3K inhibitor"),]

multi3 <- bind_rows(multi3,multi3)
multi3$moa_v1 <- c("mtor_inhibitor","pi3k_inhibitor")
labels <- labels[-which(labels$moa == "DNA dependent protein kinase inhibitor|mTOR inhibitor|PI3K inhibitor"),]

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

multi4 <- labels[which(labels$moa == "mTOR inhibitor|PI3K inhibitor"),]
multi4 <- bind_rows(multi4,multi4)
multi4$moa_v1 <- c("mtor_inhibitor","pi3k_inhibitor")
labels <- labels[-which(labels$moa == "mTOR inhibitor|PI3K inhibitor"),]
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# pi3k

pi3k_inhibitor <- c("PI3K inhibitor","DNA protein kinase inhibitor|PI3K inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% pi3k_inhibitor)] <- "pi3k_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

#plk

plk_inhibitor <- c("PLK inhibitor","cell cycle inhibitor|PLK inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% plk_inhibitor)] <- "plk_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

multiples <- bind_rows(multi,multi2,multi3,multi4)

saveRDS(multiples,"data/cmap/labels_multiples_5_5.rds")
saveRDS(labels,"data/cmap/labels_5_5.rds")
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))
multiples <- readRDS("data/cmap/labels_multiples_5_5.rds")
# cdk inhibitor

cdk <- c("CDK inhibitor","CDK inhibitor|cell cycle inhibitor|MCL1 inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% cdk)] <- "cdk_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# cdk + tyrosine

id5 <- which(labels$moa == "CDK inhibitor|tyrosine kinase inhibitor")
multi5 <- labels[id5,]
multi5 <- bind_rows(multi5,multi5)
multi5$moa_v1 <- c("cdk_inhibitor","tyrosine_kinase_inhibitor")

labels <- labels[-id5,]
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))
multiples <- bind_rows(multiples,multi5)


# cdk + glucogen synthase

id6 <- which(labels$moa == "CDK inhibitor|glycogen synthase kinase inhibitor")
multi6 <- labels[id6,]
multi6 <- bind_rows(multi6,multi6)
multi6$moa_v1 <- c("cdk_inhibitor","glycogen_synthase_kinase_inhibitor")

labels <- labels[-id6,]
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))
multiples <- bind_rows(multiples,multi6)


saveRDS(multiples,"data/cmap/labels_multiples_5_5.rds")
saveRDS(labels,"data/cmap/labels_5_5.rds")

labels <- readRDS("data/cmap/labels_5_5.rds")
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))
multiples <- readRDS("data/cmap/labels_multiples_5_5.rds")

# mapk and mek inhibitor

id <- which(labels$moa == "MAP kinase inhibitor|MEK inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("mek_inhibitor","mapk_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "MEK inhibitor")] <- "mek_inhibitor"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

#mapk inhibitors

mapk_inhibitor <- c("p38 MAPK inhibitor","MAP kinase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% mapk_inhibitor)] <- "mapk_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# nfkb and akt inhibitors

id <- which(labels$moa == "AKT inhibitor|differentiation inducer|NFKB pathway inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("akt_inhibitor","nfkb_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "NFkB pathway inhibitor")] <- "nfkb_inhibitor"
labels$moa_v1[which(labels$moa == "AKT inhibitor")] <- "akt_inhibitor"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# nfkb and ikk inhibitor

id <- which(labels$moa == "IKK inhibitor|NFkB pathway inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("ikk_inhibitor","nfkb_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "IKK inhibitor")] <- "ikk_inhibitor"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# nfkb and proteasome inhibitor

id <- which(labels$moa == "NFkB pathway inhibitor|proteasome inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("proteasome_inhibitor","nfkb_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "proteasome inhibitor")] <- "proteasome_inhibitor"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

#nfkb and tp 53 activator

id <- which(labels$moa == "cytokine production inhibitor|NFkB pathway inhibitor|TP53 activator")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("tp53_activator","nfkb_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "proteasome inhibitor")] <- "proteasome_inhibitor"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

#protein synthesis inhibitor and p450 inhibitor

id <- which(labels$moa == "cytochrome P450 inhibitor|protein synthesis inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("p450_inhibitor","protein_synthesis_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "protein synthesis inhibitor")] <- "protein_synthesis_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# retinoid receptor agonist_ligand

retinoid <- c("retinoid receptor agonist","apoptosis stimulant|retinoid receptor agonist","retinoid receptor ligand","retinoid receptor binder")
labels$moa_v1[which(labels$moa_v1 %in% retinoid)] <- "retinoid_agonist_ligand"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))


# microtubule and tubulin polymerization inhibitors

microtubule_inhibitor <- c("tubulin polymerization inhibitor",
                           "microtubule inhibitor|tubulin polymerization inhibitor",
                           "microtubule inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% microtubule_inhibitor)] <- "microtubule_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

#bacterial dna inhibitor

bact_dna <- c("bacterial DNA gyrase inhibitor","bacterial DNA inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% bact_dna)] <- "bacterial_dna_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# dopamine agonist and serotonin antagonist

id <- which(labels$moa == "dopamine receptor agonist|serotonin receptor antagonist")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("dopamine_receptor_agonist","serotonin_receptor_antagonist")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "dopamine receptor agonist")] <- "dopamine_receptor_agonist"
labels$moa_v1[which(labels$moa == "serotonin receptor antagonist")] <- "serotonin_receptor_antagonist"

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count))

# serotonin agonist 
sero_agonist <- c("serotonin receptor agonist","serotonin receptor partial agonist")
labels$moa_v1[which(labels$moa_v1 %in% sero_agonist)] <- "serotonin_receptor_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#mTOR inhibitor
labels$moa_v1[which(labels$moa == "mTOR inhibitor")] <- "mtor_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# opioid receptor agonist and antagonist

id <- which(labels$moa == "opioid receptor agonist|opioid receptor antagonist")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("opioid_receptor_agonist","opioid_receptor_antagonist")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "opioid receptor antagonist")] <- "opioid_receptor_antagonist"
labels$moa_v1[which(labels$moa == "opioid receptor agonist")] <- "opioid_receptor_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#phosphodiesterase and leukotriene inhibitors

id <- which(labels$moa == "leukotriene receptor antagonist|phosphodiesterase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("leukotriene_receptor_antagonist","phosphodiesterase_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "phosphodiesterase inhibitor")] <- "phosphodiesterase_inhibitor"
labels$moa_v1[which(labels$moa == "leukotriene receptor antagonist")] <- "leukotriene_receptor_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# glycogen synthesis and lipoxygenase

id <- which(labels$moa == "glycogen synthase kinase inhibitor|lipoxygenase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("glycogen_synthase_kinase_inhibitor","lipoxygenase_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "glycogen synthase kinase inhibitor")] <- "glycogen_synthase_kinase_inhibitor"
labels$moa_v1[which(labels$moa == "lipoxygenase inhibitor")] <- "lipoxygenase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#leukotriene and lipoxygenase

id <- which(labels$moa == "leukotriene receptor antagonist|lipoxygenase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("leukotriene_receptor_antagonist","lipoxygenase_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]

moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#phosphodiesterase inhibitor

phosphodiesterase <- c("phosphodiesterase_inhibitor","phosphodiesterase inhibitor|platelet aggregation inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% phosphodiesterase)] <- "phosphodiesterase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# calcium channel blocker

calcium <- c("calcium channel blocker","T-type calcium channel blocker")
labels$moa_v1[which(labels$moa_v1 %in% calcium)] <- "calcium_channel_blocker"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# leucine raf inhibitors

id <- which(labels$moa == "leucine rich repeat kinase inhibitor|RAF inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("leucine_inhibitor","leucine_inhibitor","raf_inhibitor","raf_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "leucine rich repeat kinase inhibitor")] <- "leucine_inhibitor"
labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


saveRDS(multiples,"data/cmap/labels_multiples_8_5.rds")
saveRDS(labels,"data/cmap/labels_8_5.rds")

multiples <- readRDS("data/cmap/labels_multiples_8_5.rds")
labels <- readRDS("data/cmap/labels_8_5.rds")
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# raf + flt3 + pdgfr + kit + vegfr

raf <- "raf_inhibitor"
flt <- "flt_inhibitor"
pdgfr <- "pdgfr_inhibitor"
vegfr <- "vegfr_inhibitor"
kit <- "kit_inhibitor"
fgfr <- "fgfr_inhibitor"
src <- "src_inhibitor"
bcr <- "bcr_abl_inhibitor"
jak <- "jak_inhibitor"
egfr <- "egfr_inhibitor"
gfr <- "gfr_inhibitor"
trk <- "tyrosine_kinase_inhibitor"
cmet <- "cmet_inhibitor"
pkc <- "pkc_inhibitor"

id <- which(labels$moa == "FLT3 inhibitor|KIT inhibitor|PDGFR tyrosine kinase receptor inhibitor|RAF inhibitor|RET tyrosine kinase inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi,multi,multi)
multi$moa_v1 <- c(flt,kit,pdgfr,raf,vegfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "leucine rich repeat kinase inhibitor")] <- "leucine_inhibitor"
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt+bcr+jak

id <- which(labels$moa == "Aurora kinase inhibitor|Bcr-Abl kinase inhibitor|FLT3 inhibitor|JAK inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(flt,bcr,jak)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt + bcr + pdgfr

id <- which(labels$moa == "Bcr-Abl kinase inhibitor|FLT3 inhibitor|PDGFR tyrosine kinase receptor inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(flt,bcr,pdgfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#flt + egfr + fgfr + pdgfr + vegfr

id <- which(labels$moa == "EGFR inhibitor|FGFR inhibitor|FLT3 inhibitor|PDGFR tyrosine kinase receptor inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi,multi,multi)
multi$moa_v1 <- c(flt,egfr,pdgfr,vegfr,fgfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# cmet + flt + trk

id <- which(labels$moa == "AXL kinase inhibitor|c-Met inhibitor|FLT3 inhibitor|hepatocyte growth factor receptor inhibitor|macrophage migration inhibiting factor inhibitor|tyrosine kinase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(flt,trk,cmet)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt + gfr + jak

id <- which(labels$moa == "FLT3 inhibitor|growth factor receptor inhibitor|JAK inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(flt,gfr,jak)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt + jak 

id <- which(labels$moa == "FLT3 inhibitor|JAK inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(flt,jak)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt + kit + pdgfr + vegfr

id <- which(labels$moa == "FLT3 inhibitor|KIT inhibitor|PDGFR tyrosine kinase receptor inhibitor|RET tyrosine kinase inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi,multi)
multi$moa_v1 <- c(flt,kit,pdgfr,vegfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# flt + kit + pkc

id <- which(labels$moa == "FLT3 inhibitor|KIT inhibitor|PKC inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(flt,kit,pkc)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


# pdgfr + vegfr + fgfr

id <- which(labels$moa == "FGFR inhibitor|PDGFR tyrosine kinase receptor inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(pdgfr,vegfr,fgfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#pdgfr + kit + vegfr

id <- which(labels$moa == "KIT inhibitor|PDGFR tyrosine kinase receptor inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(pdgfr,vegfr,kit)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#pdgfr bcr kit src trk

id <- which(labels$moa == "Bcr-Abl kinase inhibitor|ephrin inhibitor|KIT inhibitor|PDGFR tyrosine kinase receptor inhibitor|SRC inhibitor|tyrosine kinase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi,multi,multi)
multi$moa_v1 <- c(pdgfr,bcr,kit,src,trk)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# pdgfr bcr kit trk 
id <- which(labels$moa == "Bcr-Abl kinase inhibitor|KIT inhibitor|PDGFR tyrosine kinase receptor inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi,multi)
multi$moa_v1 <- c(pdgfr,bcr,kit,trk)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

labels$moa_v1[which(labels$moa == "PDGFR tyrosine kinase receptor inhibitor")] <- pdgfr
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

raf <- "raf_inhibitor"
flt <- "flt_inhibitor"
pdgfr <- "pdgfr_inhibitor"
vegfr <- "vegfr_inhibitor"
kit <- "kit_inhibitor"
fgfr <- "fgfr_inhibitor"
src <- "src_inhibitor"
bcr <- "bcr_abl_inhibitor"
jak <- "jak_inhibitor"
egfr <- "egfr_inhibitor"
gfr <- "gfr_inhibitor"
trk <- "tyrosine_kinase_inhibitor"
cmet <- "cmet_inhibitor"
pkc <- "pkc_inhibitor"

# bcr + src + vegfr

id <- which(labels$moa == "Abl kinase inhibitor|SRC inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi,multi)
multi$moa_v1 <- c(vegfr,bcr,src)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#vegfr fgfr

id <- which(labels$moa == "FGFR inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(vegfr,fgfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "FLT3 inhibitor")] <- flt
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#vegfr kit 

id <- which(labels$moa == "KIT inhibitor|VEGFR inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(vegfr,kit)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "VEGFR inhibitor")] <- vegfr
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


#bcr src

id <- which(labels$moa == "Abl kinase inhibitor|Bcr-Abl kinase inhibitor|SRC inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(src,bcr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "SRC inhibitor")] <- src
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#src egfr 

id <- which(labels$moa == "EGFR inhibitor|SRC inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(src,egfr)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "SRC inhibitor")] <- src
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


labels$moa_v1[which(labels$moa == "Abl kinase inhibitor|Bcr-Abl kinase inhibitor")] <- bcr
labels$moa_v1[which(labels$moa == "Bcr-Abl kinase inhibitor")] <- bcr
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


#egfr jak

id <- which(labels$moa == "EGFR inhibitor|JAK inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(egfr,jak)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "JAK inhibitor")] <- jak
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

#egfr pkc

id <- which(labels$moa == "EGFR inhibitor|PKC inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c(egfr,pkc)
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "EGFR_inhibitor")] <- egfr
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


labels$moa_v1[which(labels$moa == "growth factor receptor inhibitor")] <- gfr
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

trk_in <- c("tyrosine kinase inhibitor","protein tyrosine kinase inhibitor","protein tyrosine kinase inhibitor|tyrosine kinase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% trk_in)] <- trk
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


labels$moa_v1[which(labels$moa == "c-Met inhibitor")] <- cmet
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


labels$moa_v1[which(labels$moa == "PKC inhibitor")] <- pkc
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# estrogen antagonist

estro_antagonist <- c("estrogen receptor antagonist","estrogen receptor antagonist|selective estrogen receptor modulator (SERM)")
labels$moa_v1[which(labels$moa_v1 %in% estro_antagonist)] <- "estrogen_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# estrogen antagonist and progesterone agonist

id <- which(labels$moa == "estrogen receptor antagonist|progesterone receptor agonist")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("estrogen_antagonist","progesterone_agonist")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
#labels$moa_v1[which(labels$moa == "EGFR_inhibitor")] <- egfr
#labels$moa_v1[which(labels$moa == "RAF inhibitor")] <- "raf_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

labels$moa_v1[which(labels$moa == "estrogen receptor agonist")] <- "estrogen_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


ribo_nuc <- c("ribonucleotide reductase inhibitor","adenosine deaminase inhibitor|ribonucleotide reductase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% ribo_nuc)] <- "ribonucleotide_reductase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


p450 <- c("cytochrome P450 inhibitor","cytochrome P450 inhibitor|imidazoline receptor ligand")
labels$moa_v1[which(labels$moa_v1 %in% p450)] <- "p450_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

pota_activ <- c("potassium channel activator","nitric oxide donor|potassium channel activator")
labels$moa_v1[which(labels$moa_v1 %in% pota_activ)] <- "potassium_activator"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


labels$moa_v1[which(labels$moa == "progesterone receptor agonist")] <- "progesterone_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


bact_ribo <- c("bacterial 50S ribosomal subunit inhibitor","bacterial 30S ribosomal subunit inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% bact_ribo)]<- "bacterial_ribosomal_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

antibacterial <- c("bacterial_ribosomal_inhibitor","antibacterial agent","bacterial antifolate")
labels$moa_v1[which(labels$moa_v1 %in% antibacterial)]<- "antibacterial"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# dna methyltransferase and trpv agonist

id <- which(labels$moa == "aldehyde dehydrogenase inhibitor|DNA methyltransferase inhibitor|TRPV agonist")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("dna_methyltransferase_inhibitor","trpv_agonist")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "DNA methyltransferase inhibitor")] <- "dna_methyltransferase_inhibitor"
labels$moa_v1[which(labels$moa == "TRPV agonist")] <- "trpv_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


ppar <- c("PPAR receptor agonist","insulin sensitizer|PPAR receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% ppar)]<- "ppar_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


prosta_agonist <- c("prostanoid receptor agonist","platelet aggregation inhibitor|prostanoid receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% prosta_agonist)]<- "prostanoid_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


syk <- c("SYK inhibitor","mediator release inhibitor|SYK inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% syk)]<- "syk_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

bcl <- c("BCL inhibitor","BCL inhibitor|MCL1 inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% bcl)]<- "bcl_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

chk <- c("CHK inhibitor","CHK inhibitor|cytochrome P450 activator|indoleamine 2,3-dioxygenase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% chk)]<- "chk_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


# dehydrogenase inositol

id <- which(labels$moa == "dehydrogenase inhibitor|inositol monophosphatase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("dehydrogenase_inhibitor","dehydrogenase_inhibitor","dehydrogenase_inhibitor","inositol_inhibitor","inositol_inhibitor","inositol_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "inositol monophosphatase inhibitor")] <- "inositol_inhibitor"

labels$moa_v1[which(labels$moa == "11-beta hydroxysteroid dehydrogenase inhibitor|FXR agonist")] <- "dehydrogenase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


vitamin_d <- c("vitamin D receptor agonist","FXR antagonist|vitamin D receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% vitamin_d)]<- "vitamin_d_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

adenosine <- c("adenosine receptor agonist|cAMP stimulant","adenosine receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% adenosine)]<- "adenosine_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

anticancer <- c("anticancer agent","anticancer agent|aryl hydrocarbon receptor antagonist")
labels$moa_v1[which(labels$moa_v1 %in% anticancer)]<- "anticancer"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

faah <- c("FAAH inhibitor","FAAH inhibitor|FAAH reuptake inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% faah)]<- "faah_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

antifungal <- c("fungal squalene epoxidase inhibitor","other antifungal")
labels$moa_v1[which(labels$moa_v1 %in% antifungal)]<- "antifungal"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


snri <- c("norepinephrine reputake inhibitor","adrenergic inhibitor|norepinephrine reuptake inhibitor|serotoninâ€“norepinephrine reuptake inhibitor (SNRI)",
          "norepinephrine inhibitor|norepinephrine reuptake inhibitor|serotonin receptor antagonist|serotoninâ€“norepinephrine reuptake inhibitor (SNRI)",
          "norepinephrine reuptake inhibitor","serotoninâ€“norepinephrine reuptake inhibitor (SNRI)")
labels$moa_v1[which(labels$moa_v1 %in% snri)]<- "snri"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


transported <- c("sodium/potassium/chloride transporter inhibitor","sodium/glucose cotransporter inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% transported)]<- "sodium_potassium_chloride_transporter_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# acetylcholine agonist

acetyl <- c("acetylcholine release enhancer","acetylcholine receptor agonist")
labels$moa_v1[which(labels$moa_v1 %in% acetyl)]<- "acetylcholine_agonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# androgen

androgen <- c("androgen receptor antagonist","androgen receptor modulator")
labels$moa_v1[which(labels$moa_v1 %in% androgen)]<- "androgen_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

# arrhythmia

arrhy <- c("antiarrhythmic","antiarrhythmic medication")
labels$moa_v1[which(labels$moa_v1 %in% arrhy)]<- "anti_arrhythmic"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


saveRDS(multiples,"data/cmap/labels_multiples_10_5.rds")
saveRDS(labels,"data/cmap/labels_10_5.rds")

library(tidyverse)
multiples <- readRDS("data/cmap/labels_multiples_10_5.rds")
labels <- readRDS("data/cmap/labels_10_5.rds")
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))


anticancer <- c("antineoplastic agent","anticancer")
labels$moa_v1[which(labels$moa_v1 %in% anticancer)]<- "anticancer"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

trpv_ant <- c("aromatase inhibitor|TRPV antagonist","TRPV antagonist")
labels$moa_v1[which(labels$moa_v1 %in% trpv_ant)]<- "trpv_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

moa_copy <- moa_v1

#benzo ago and HDAC inhibitor

id <- which(labels$moa == "benzodiazepine receptor agonist|HDAC inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("hdac_inhibitor","benzo_agonist")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "benzodiazepine receptor agonist")] <- "benzo_agonist"

labels$moa_v1[which(labels$moa == "HDAC inhibitor")] <- "hdac_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


choli <- c("acetylcholinesterase inhibitor","cholinesterase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% choli)]<- "acetylcholinesterase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


kin_inh <- c("protein kinase inhibitor","diacylglycerol kinase inhibitor|protein kinase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% kin_inh)]<- "kinase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1

peptidase <- c("dehydropeptidase inhibitor","dipeptidyl peptidase inhibitor","peptidase inhibitor")
labels$moa_v1[which(labels$moa_v1 %in% peptidase)]<- "peptidase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1

#dna pol and thymid

id <- which(labels$moa == "DNA directed DNA polymerase inhibitor|thymidylate synthase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("dna_pol_inhibitor","thymidylate_synthase_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "thymidylate synthase inhibitor")] <- "thymidylate_synthase_inhibitor"

labels$moa_v1[which(labels$moa == "DNA polymerase inhibitor")] <- "dna_pol_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1

gaba <- c("GABA receptor antagonist","GABA receptor antagonist|GABA receptor modulator",
          "GABA receptor antagonist|TP53 activator","GABA receptor modulator","GABA receptor negative allosteric modulator")

labels$moa_v1[which(labels$moa_v1 %in% gaba)]<- "gaba_antagonist"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1

# hiv integrase + nucleoside reverse transcriptase

id <- which(labels$moa == "HIV integrase inhibitor|nucleoside reverse transcriptase inhibitor")
multi <- labels[id,]
multi <- bind_rows(multi,multi)
multi$moa_v1 <- c("hiv_integrase","nucleoside_reverse_transcriptase_inhibitor")
multiples <- bind_rows(multiples,multi)

labels <- labels[-id,]
labels$moa_v1[which(labels$moa == "nucleoside reverse transcriptase inhibitor")] <- "nucleoside_reverse_transcriptase_inhibitor"

labels$moa_v1[which(labels$moa == "HIV integrase inhibitor")] <- "hiv_integrase"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


labels$moa_v1[which(labels$moa == "leucine rich repeat kinase inhibitor")] <- "leucine_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


muscle <- c("muscle relaxant","myorelaxant")

labels$moa_v1[which(labels$moa_v1 %in% muscle)]<- "muscle_relaxant"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1

phosphatase <- c("mitochondrial DNA polymerase inhibitor|phosphatase inhibitor","protein phosphatase inhibitor")

labels$moa_v1[which(labels$moa_v1 %in% phosphatase)]<- "phosphatase_inhibitor"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


rna_dna <- c("dna_rna_damage","RNA polymerase inhibitor")

labels$moa_v1[which(labels$moa_v1 %in% rna_dna)]<- "dna_rna_damage"
moa_v1 <- labels %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
moa_copy <- moa_v1


saveRDS(multiples,"data/cmap/labels_multiples_10_5.rds")
saveRDS(labels,"data/cmap/labels_10_5.rds")
multiples <- readRDS("data/cmap/labels_multiples_10_5.rds")
labels <- readRDS("data/cmap/labels_10_5.rds")
all <- bind_rows(labels,multiples)
all_moa_v1 <- all %>% group_by(moa_v1) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))

saveRDS(all,"data/cmap/labels_first_pass.rds")
all_orig <- readRDS("data/cmap/labels.rds")

all_orig <- all_orig %>% mutate(target = str_split(string = target,pattern = "\\|")) %>% unnest()

targets <- all_orig %>% group_by(target) %>% summarise(count = n()) %>% arrange(desc(count)) %>% mutate(cs = cumsum(count))
