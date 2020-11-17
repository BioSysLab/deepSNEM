library(tidyverse)
library(topGO)
# get all unique gene symbols from all the graphs unweighted

# file info
file_info <- readRDS("data/graph_info_df/file_info_nodups.rds")

# dup file

file_info_dups <- readRDS("data/graph_info_df/file_info_dups.rds")

proteins <- NULL

for (i in 1:nrow(file_info)) {
  graph <- read.csv(file_info$files_combined[i])
  symbols <- unique(c(as.character(graph$node1),as.character(graph$node2)))
  proteins <- c(proteins,symbols)
  proteins <- unique(proteins)
}


for (i in 1:nrow(file_info_dups)) {
  graph <- read.csv(file_info_dups$files_combined[i])
  symbols <- unique(c(as.character(graph$node1),as.character(graph$node2)))
  proteins <- c(proteins,symbols)
  proteins <- unique(proteins)
}

proteins <- as.data.frame(proteins)

proteins <- proteins %>% mutate(proteins = str_replace_all(string = proteins,pattern = "~",replacement = "-"))

genes <- factor(x = rep(1,nrow(proteins)),levels = c(0,1))
names(genes) <- as.character(proteins$proteins)

GOobject_BP <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                   ID = "SYMBOL", nodeSize = 100)

term.genes_bp <- genesInTerm(GOobject_BP, GOobject_BP@graph@nodes)

GOobject_CC <- new("topGOdata",ontology = "CC", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                   ID = "SYMBOL", nodeSize = 100)

term.genes_cc <- genesInTerm(GOobject_CC, GOobject_CC@graph@nodes)


GOobject_MF <- new("topGOdata",ontology = "MF", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                   ID = "SYMBOL", nodeSize = 100)

term.genes_mf <- genesInTerm(GOobject_MF, GOobject_MF@graph@nodes)

bp <- matrix(0,nrow = nrow(proteins), ncol = length(term.genes_bp))
for (i in 1:nrow(proteins)) {
  for (j in 1:length(term.genes_bp)) {
    if (any(proteins$proteins[i] %in% term.genes_bp[[j]])) {
      bp[i,j] <- 1
    }
  }
}
x <- apply(X = bp,MARGIN = 1,sum)
x <- as.data.frame(x)

mf <- matrix(0,nrow = nrow(proteins), ncol = length(term.genes_mf))
for (i in 1:nrow(proteins)) {
  for (j in 1:length(term.genes_mf)) {
    if (any(proteins$proteins[i] %in% term.genes_mf[[j]])) {
      mf[i,j] <- 1
    }
  }
}
x_mf <- apply(X = mf,MARGIN = 1,sum)
x_mf <- as.data.frame(x_mf)

cc <- matrix(0,nrow = nrow(proteins), ncol = length(term.genes_cc))
for (i in 1:nrow(proteins)) {
  for (j in 1:length(term.genes_cc)) {
    if (any(proteins$proteins[i] %in% term.genes_cc[[j]])) {
      cc[i,j] <- 1
    }
  }
}

x_cc <- apply(X = cc,MARGIN = 1,sum)
x_cc <- as.data.frame(x_cc)

proteins <- proteins %>% mutate(proteins = str_replace_all(string = proteins,pattern = "-",replacement = "~"))
rownames(bp) <- as.character(proteins$proteins)
rownames(cc) <- as.character(proteins$proteins)
rownames(mf) <- as.character(proteins$proteins)

colnames(bp) <- as.character(GOobject_BP@graph@nodes)
colnames(cc) <- as.character(GOobject_CC@graph@nodes)
colnames(mf) <- as.character(GOobject_MF@graph@nodes)

write.csv(bp,"data/prot_embeddings/unweighted/bp_features.csv",row.names = T)
write.csv(mf,"data/prot_embeddings/unweighted/mf_features.csv",row.names = T)
write.csv(cc,"data/prot_embeddings/unweighted/cc_features.csv",row.names = T)
write.csv(proteins,"data/prot_embeddings/unweighted/proteins.csv")


# weighted protein features

