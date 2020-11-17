drug_sigs_per_line <- function(cell_line,sig_info,sig_metrics) {
  
  # This functions takes as input a cell line, the CMAP signature info and metrics
  # and returns a dataframe with 1 signature per unique compound.
  # Multiple signatures of the same drug are filtered by taking into account the quality of the signatures.
  # To calculate the quality a combination of is_exemplar, TAS and number of replicates is used.
  # cell_line character of cell line
  # sig_info dataframe of GSE info
  # sig_metrics dataframe of GSE metrics
  
  library(tidyverse)
  options(warn =- 1)
  
  
  cell <- sig_info %>%
    filter(cell_id == cell_line) %>%
    filter(pert_type == "trt_cp") %>%
    group_by(pert_iname) %>%
    mutate(count = n_distinct(sig_id)) %>%
    ungroup()
  
  print(paste0('the unique drugs for ',cell_line,' are ',length(unique(cell$pert_iname))))
  
  ## add the signature metrics
  
  cell <- left_join(cell,sig_metrics)
  
  ## keep the drugs that we have only 1 signature for this cell line
  
  cell_singles <- cell %>%
    filter(count == 1) %>%
    dplyr::select(-count)
  
  print(paste0('the drugs that have only 1 signature for ',cell_line,' are ',length(unique(cell_singles$pert_iname))))
  
  cell_singles$pert_itime <- factor(cell_singles$pert_itime)
  print("time summary")
  print(summary(cell_singles$pert_itime))
  cell_singles$pert_idose <- factor(cell_singles$pert_idose)
  print("dose summary")
  print(summary(cell_singles$pert_idose))
  
  ## add quality column to single perturbations
  
  cell_singles$quality <- 100
  
  cell_singles <- cell_singles %>%
    mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
           quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
           quality = factor(quality))
  
  print("summary of the quality of drugs with only 1 signature")
  print(summary(cell_singles$quality))
  
  ## keep the multiple signature drugs in cell
  
  cell<- anti_join(cell,cell_singles)
  
  ### add priorities to the multiple signatures
  
  cell$priority <- 100
  cell <- cell %>%
    mutate(priority = if_else(pert_dose == "10.0" & pert_time == 24,true = 1,false = priority),
           priority = if_else(pert_idose == "5 ÂµM" & pert_time == 24,true = 2,false = priority),
           priority = if_else(pert_idose != "5 ÂµM" & pert_dose != "10.0" & pert_time == 24,true = 3,false = priority),
           priority = if_else(pert_dose == "10.0" & pert_time == 6,true = 4,false = priority),
           priority = if_else(pert_idose == "5 ÂµM" & pert_time == 6,true = 5,false = priority),
           priority = if_else(pert_idose != "5 ÂµM" & pert_dose != "10.0" & pert_time == 6,true = 6,false = priority),
           priority = factor(priority))
  
  print("priorities for drugs with multiple signatures")
  print(summary(cell$priority))
  ### add quality to the multiple signatures
  
  cell$quality <- 100
  
  cell <- cell %>%
    mutate(quality = if_else(is_exemplar == 1 & tas > 0.4 & distil_nsample>=2 ,true = 1,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample>2 ,true = 2,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.2 & tas<=0.4 & distil_nsample <=2 ,true = 3,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample>2 ,true = 4,false = quality),
           quality = if_else(is_exemplar == 1 & tas > 0.1 & tas<=0.2 & distil_nsample <= 2 ,true = 5,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample > 2 ,true = 6,false = quality),
           quality = if_else(is_exemplar == 1 & tas < 0.1 & distil_nsample <= 2 ,true = 7,false = quality),
           quality = if_else(is_exemplar == 0 ,true = 8,false = quality),
           quality = factor(quality))
  
  
  print("summary of the quality of drugs with multiple signatures")
  print(summary(cell$quality))
  
  
  print(paste0('the drugs that have Multiple signatures for ',cell_line,' are ',length(unique(cell$pert_iname))))
  
  
  
  
  #### clean them based on quality for each drug and then solve the equalities with max tas
  
  
  
  
  cell_cleaned <- cell %>%
    group_by(pert_iname) %>%
    filter(quality == min(as.numeric(quality))) %>%
    filter(tas == max(tas)) %>%
    ungroup %>%
    dplyr::select(-c(count,priority))
  
  
  
  cell_final <- bind_rows(cell_cleaned,cell_singles)
  
  print("summary of final quality of signatures")
  print(summary(cell_final$quality))
  
  return(cell_final)
}
