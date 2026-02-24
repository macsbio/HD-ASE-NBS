#Clear working space
rm(list = ls(all.names = T))

# Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)}
if(!require(ggsci)){
  install.packages("ggsci")
  library(ggsci)}
if(!require(qvalue)){
  install.packages("qvalue")
  library(qvalue)}
if(!require(dplyr)) {
  install("dplyr")
  library(dplyr) }
if(!require(tidyr)) {
  install("tidyr")
  library(tidyr) }
if(!require(rJava)) {
  install.packages("rJava")
  library(rJava) }
if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)}
if(!require(readxl)){
  install.packages("readxl")
  library(readxl)}
if(!require(xlsx)){
  install.packages("xlsx")
  library(xlsx)}
if(!require(pROC)){
  install.packages("pROC")
  library(pROC)}
if(!require(BiocManager)) {
  install.packages("BiocManager") }
if(!require(biomaRt)){
  BiocManager::install("biomaRt")
  library(biomaRt)}
if(!require(limma)) {
  BiocManager::install("limma")
  library(limma) }
if(!require(edgeR)) {
  BiocManager::install("edgeR") }
library(edgeR)
if(!require(RColorBrewer)) {
  BiocManager::install("RColorBrewer") }
library(RColorBrewer)
if(!require(ggrepel)) {
  BiocManager::install("ggrepel") }
library(ggrepel)
if(!require(pcaMethods)) {
  BiocManager::install("pcaMethods") }
library(pcaMethods)
if(!require(org.Hs.eg.db)){
  BiocManager::install("org.Hs.eg.db") }
library(org.Hs.eg.db)
if(!require(topGO)){
  BiocManager::install("topGO") }
library(topGO)
if(!require("RCy3")) {
  BiocManager::install("RCy3") }
library(RCy3)

# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
if(!dir.exists(paste0(DATA.DIR, "/Output"))) {
  dir.create(paste0(DATA.DIR, "/Output"))
} 
RES.DIR <- paste0(DATA.DIR, "/Output") 

options(stringsAsFactors = F)

# Load processed ASEReadCounter output data
setwd(DATA.DIR)


idnetwork <- read.table("HTT_900_cc_matrix.csv", header = T, sep = ",")
rownames(idnetwork) <- idnetwork$X
idnetwork <- idnetwork[,-1]
all(rownames(idnetwork) == colnames(idnetwork))

idnetwork <- as.matrix(idnetwork)

cytoscapePing()

idgraph <- graph_from_adjacency_matrix(idnetwork,mode = "undirected", weighted = "weight")
createNetworkFromIgraph(idgraph, "patient", "patient")
idcluster <- read.table("HTT_900_wasp_cluster_assignments.csv", header = F, sep = ",")
colnames(idcluster) <- c("id","cluster")
rownames(idcluster) <- idcluster$id
loadTableData(data = idcluster, network = "patient")
deleteSelfLoops(network = "patient")


#save.image("PYNBS_NETWORK.RData")
# load("PYNBS_NETWORK.RData")
# 
# prop07 <- read.table("ResultsDCM_ASE_prop_kernel.csv", header = T, sep = ",")
# rownames(prop07) <- prop07$X
# prop07 <- prop07[,-1]
# 
# prop07 <- t(prop07)
# prop07 <- as.data.frame(prop07)
# prop07_95 <- prop07[rowSums(prop07)>=quantile(rowSums(prop07), prob = 0.05),]
# 
# percImbGene <- as.data.frame(prop07)
# percImbGene$perc1 <- NA
# percImbGene$perc2 <- NA
# percImbGene$perc3 <- NA
# percImbGene$perc4 <- NA
# percImbGene$sum1 <- NA
# percImbGene$sum2 <- NA
# percImbGene$sum3 <- NA
# percImbGene$sum4 <- NA
# percImbGene$prop1 <- NA
# percImbGene$prop2 <- NA
# percImbGene$prop3 <- NA
# percImbGene$prop4 <- NA
# percImbGene <- percImbGene[,c("perc1", "perc2", "perc3", "perc4",
#                               "sum1","sum2","sum3","sum4")]
# 
# clusters <- read.table("cluster_assignment.txt", sep = "\t")
# colnames(clusters) <- c("SID","Cluster")
# 
# networkdataT <- as.data.frame(t(networkdata))
# 
# clusters$perc[clusters$Cluster==1] <- "perc1"
# clusters$perc[clusters$Cluster==2] <- "perc2"
# clusters$perc[clusters$Cluster==3] <- "perc3"
# clusters$perc[clusters$Cluster==4] <- "perc4"
# rownames(clusters) <- clusters$SID
# 
# for (i in unique(clusters$Cluster)) {
#   binary <- networkdataT[rownames(networkdataT) %in% rownames(prop07_95),colnames(networkdataT) %in% 
#                            clusters$SID[clusters$Cluster==i]]
#   for (g in rownames(percImbGene)) {
#     percImbGene[g,paste0("perc",i)] <- (rowSums(binary[g,])/ncol(binary))
#     percImbGene[g,paste0("sum",i)] <- rowSums(binary[g,])
#   }
# }
# 
# for (i in unique(clusters$Cluster)) {
#   propagation <- prop07_95[rownames(networkdataT) %in% rownames(prop07_95),colnames(networkdataT) %in% 
#                              clusters$SID[clusters$Cluster==i]]
#   for (g in rownames(percImbGene)) {
#     percImbGene[g,paste0("prop",i)] <- rowSums(propagation[g,])
#   } 
# }
# 
# percImbGene$ensembl_gene_id <- rownames(percImbGene)
# write.table(percImbGene, "imbalance_%_cluster.txt", col.names = T, row.names = T, quote = F, sep = "\t")
# 
# save.image("PYNBS_NETWORK.RData")
# 
# # MANN-WHITNEY U TEST PER CLUSTER VS. REST HEAT PROPAGATION 
# mwVals <- list()
# mwConf <- list()
# for (i in unique(clusters[!is.na(clusters$Cluster),"Cluster"])) {
#   newLabels <- clusters
#   newLabels <- na.omit(newLabels)
#   newLabels$Cluster <- as.numeric(newLabels$Cluster)
#   newLabels[!newLabels$Cluster==i, "Cluster"] <- 0
#   newLabels[newLabels$Cluster==i, "Cluster"] <- 1
#   newLabels$Cluster <- as.factor(newLabels$Cluster)
#   mwVals[[paste0("Cluster ", i)]] <- apply(prop07_95, 1, function(x) {
#     testData <- split(as.vector(as.matrix(x)), newLabels[, "Cluster"])
#     mw <- NA
#     mw <- wilcox.test(testData[[1]], 
#                       testData[[2]], na.rm = T)$p.value
#     return(mw)
#   })
#   mwConf[[paste0("Cluster ", i)]] <- apply(prop07_95, 1, function(x) {
#     testData <- split(as.vector(as.matrix(x)), newLabels[, "Cluster"])
#     int <- NA
#     int <- wilcox.test(testData[[1]], 
#                        testData[[2]], na.rm = T, conf.int = T)$estimate
#     return(int)
#   })
#   names(mwVals[[paste0("Cluster ", i)]]) <- rownames(prop07_95)
#   names(mwConf[[paste0("Cluster ", i)]]) <- rownames(prop07_95)
# }
# 
# networkGenesRankMW <- as.data.frame(mwVals)
# networkGenesRankMW1 <- na.omit(networkGenesRankMW)
# networkGenesRankMW1$ensembl_gene_id <- rownames(networkGenesRankMW1)
# 
# networkGenesEff <- as.data.frame(mwConf)
# networkGenesEff <- na.omit(networkGenesEff)
# colnames(networkGenesEff) <- c("Cluster2eff", "Cluster3eff", "Cluster4eff", "Cluster1eff")
# networkGenesEff$ensembl_gene_id <- rownames(networkGenesEff)
# 
# networkGenesMW <- full_join(networkGenesEff,networkGenesRankMW1)
# 
# # Take highest ASE hit per person per gene
# library(dplyr)
# absScores1 <- absScores %>% distinct(SID, ensembl_gene_id, .keep_all= TRUE)
# absScores1 <- na.omit(absScores1)
# 
# clinData$SID <- as.character(clinData$SID)
# absScores1Cluster <- full_join(absScores1, clinData[,c("SID","Group")])
# absScores1Cluster <- na.omit(absScores1Cluster)
# 
# medianScore <- data.frame(matrix(nrow = length(unique(absScores1Cluster$ensembl_gene_id)),
#                                  ncol = 4))
# rownames(medianScore) <- unique(absScores1Cluster$ensembl_gene_id)
# for (i in 1:4) {
#   for (j in rownames(medianScore)) {
#     medianScore[j,i] <- median(na.omit(absScores1Cluster[absScores1Cluster$Group==i & 
#                                                            absScores1Cluster$ensembl_gene_id==j,])$ASE)
#   }
# }

