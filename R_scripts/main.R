library(tidyverse)
library(cmapR)
library(rhdf5)
library(CARNIVAL)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(viper)


# cmap path to lvl5
ds_path <- "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

# read cmap signature info dataframe

sig_info <- read.delim(file = "C:/Users/user/Documents/phd/GSE92742_Broad_LINCS_sig_info.txt")

# read cmap signature metrics

sig_metrics <- read.delim(file = "cmap/GSE92742_Broad_LINCS_sig_metrics.txt")

### read landmark genes

landmark <- read_tsv(file = "cmap/cmap_landmark_genes.txt")



A375 <- drug_sigs_per_line(cell_line = "A375",sig_info = sig_info,sig_metrics = sig_metrics)

cmap_A375 <- get_cmap_signatures(cmap_path_to_gctx = ds_path,sig_ids = A375$sig_id,landmark = TRUE,landmark_df = landmark)

#dorothea

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL")) # loading the viper regulons
TF_cmap_A375 <-runDoRothEA(cmap_A375, regulon=viper_regulon, confidence_level=c('A','B','C')) # Estimating TF activities


