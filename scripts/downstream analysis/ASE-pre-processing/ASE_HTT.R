#=============================================================================#
# DCM_MUMC_ASE_analysis.R                                                     #
#                                                                             #
# Version: 4.0                                                                #
# Date:  Fall 2020 - Spring 2021                                              #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# Additional: Daan van Beek; MaCSBio                                          #
# History:                                                                    #
#  1.0: Creation                                                              #
#  2.0: Updated for new data and objects                                      #
#  2.1: Removed all ideas, thoughts, etc. from start of script and merged     #
#       into notes (Onenote)                                                  #
#                                                                             #
#  3.0: Likely homozygote filtering using WES data                            #
#  3.1: Updated and streamlined SNP <-> gene functionality (fast now!)        #
#  3.2: Bug fixing in SNP <-> gene lookup functions, added LD lookup function #
#  3.3: Updated with semi-complete data (n = 80) and phen. clustering         #
#  4.0: Update with n = 87 data, absolute deviation, new plotting             #
#                                                                             #
#=============================================================================#

#-----------------------------------------------------------------------------#
# 1. Initialization
#-----------------------------------------------------------------------------#

# Some information about necessary files and column headers
# Clinical file should contain an id header and Group (for cluster number or case/control) header
# For other general column names use:
#   EnsgID header --> ensemble_gene_id
#   HGNC header --> hgnc_symbol
#   SNP id header --> rsid
#   Stat values header --> pvalue, qvalue
# Other input files are created by running ASEReadCounter_OutputProcessing
#   on the output from the ASEReadCounter pre-processing in bash/Linux
# If interested in comparisons with known related genes
#   change the name of DCM_genes.txt
#   It should have the header ensembl_gene_id and hgnc_symbol


# Load packages
if(!require(rstudioapi, quietly = TRUE)){
  install.packages("rstudioapi")}
library(rstudioapi)

if(!require(ggplot2, quietly = TRUE)){
  install.packages("ggplot2")}
library(ggplot2)

if(!require(ggsci, quietly = TRUE)){
  install.packages("ggsci")}
library(ggsci)

if(!require(dplyr, quietly = TRUE)) {
  install.packages("dplyr") }
library(dplyr)

if(!require(tidyr, quietly = TRUE)) {
  install.packages("tidyr") }
library(tidyr)

if(!require(reshape2, quietly = TRUE)){
  install.packages("reshape2") }
library(reshape2)

if(!require(readxl, quietly = TRUE)){
  install.packages("readxl") }
library(readxl)

if(!require(writexl, quietly = TRUE)){
  install.packages("writexl") }
library(writexl)

if(!require(pROC, quietly = TRUE)){
  install.packages("pROC") }
library(pROC)

if(!require(BiocManager, quietly = TRUE)) {
  install.packages("BiocManager") }
library(BiocManager)

if(!require(biomaRt, quietly = TRUE)){
  BiocManager::install("biomaRt") }
library(biomaRt)

if(!require(qvalue, quietly = TRUE)){
  install.packages("qvalue")}
library(qvalue)

if(!require(limma, quietly = TRUE)) {
  BiocManager::install("limma") }
library(limma) 

if(!require(edgeR, quietly = TRUE)) {
  BiocManager::install("edgeR") }
library(edgeR)

if(!require(RColorBrewer, quietly = TRUE)) {
  install.packages("RColorBrewer") }
library(RColorBrewer)

if(!require(ggrepel, quietly = TRUE)) {
  install.packages("ggrepel") }
library(ggrepel)

if(!require(pcaMethods, quietly = TRUE)) {
  BiocManager::install("pcaMethods") }
library(pcaMethods)

if(!require(org.Hs.eg.db, quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db") }
library(org.Hs.eg.db)

if(!require(topGO, quietly = TRUE)){
  BiocManager::install("topGO") }
library(topGO)

if(!require("RCy3", quietly = TRUE)) {
  BiocManager::install("RCy3") }
library(RCy3)

if(!require("STRINGdb", quietly = TRUE)) {
  install("STRINGdb") }
library(STRINGdb)

if(!require("igraph", quietly = TRUE)) {
  install("igraph") }
library(igraph)

if(!require("Rgraphviz", quietly = TRUE)) {
  install("Rgraphviz") }
library(Rgraphviz)

if(!require("CePa", quietly = TRUE)) {
  install("CePa") }
library(CePa)

if(!require("RCy3", quietly = TRUE)) {
  install("RCy3") }
library(RCy3)

if(!require("plyr", quietly = TRUE)) {
  install("plyr") }
library(plyr)

if(!require("qusage", quietly = TRUE)) {
  install("qusage") }
library(qusage)


#Clear working space
rm(list = ls(all.names = T))

ASEReadCounterFile <- "ASEReadCounterOutput_Processed_HTT.RData"

# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
if(!dir.exists(paste0(DATA.DIR, "/Output"))) {
  dir.create(paste0(DATA.DIR, "/Output"))
} 
RES.DIR <- paste0(DATA.DIR, "/Output") 

options(stringsAsFactors = F)

# Load processed ASEReadCounter output data
setwd(DATA.DIR)
load(ASEReadCounterFile)

# Check if all orders of columns and rows are equal across the dataframes
all(dosageData[, 1] == altRatioData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altRatioData)) # TRUE
all(dosageData[, 1] == totalCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(totalCountData)) # TRUE
all(dosageData[, 1] == refCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(refCountData)) # TRUE
all(dosageData[, 1] == altCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altCountData)) # TRUE
all(dosageData[, 1] == refRatioData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altCountData)) # TRUE

# Restructure data
rwNms <- dosageData[, 1]
dosageData   <- dosageData[, -1]
altRatioData <- altRatioData[, -1]
totalCountData <- totalCountData[, -1]
refCountData <- refCountData[, -1]
refRatioData <- refRatioData[, -1]
altCountData <- altCountData[, -1]
rownames(altRatioData) <- rownames(dosageData) <- rownames(totalCountData) <- rwNms
rownames(refCountData) <- rownames(altCountData) <- rownames(refRatioData) <- rwNms

# Calculate the ASE score according to what is common practice:
# Ranging from 0.5 for no imbalance to 1 for full imbalance
# This is achieved by taking the alternative allele count divided by the
# total allele count (a value between 0 and 1), minus 0.5 in absolute value, + 0.5.
absData <- (abs(altRatioData - 0.5) + 0.5)

# This part is only data cleaning for OUR PROJECT
colnames(altRatioData) <- colnames(dosageData) <- colnames(totalCountData) <- colnames(absData) <- gsub("Aligned.*","", colnames(absData))

colnames(altCountData) <- colnames(refCountData) <- colnames(altRatioData)
colnames(refRatioData) <- colnames(altCountData)

#-----------------------------------------------------------------------------#
# 2. Filter the data based on likely homozygosity
#-----------------------------------------------------------------------------#

# Derive threshold for calling homozygotes from ASE ratios
# In this without genotype data
# The one we use here is based on including only loci with >= 10 reads
#                                                                 per allele
# Can be derived with other thresholds as well
#-----------------------------------------------------------------------------#
altCountData <- (altRatioData * totalCountData)
refCountData <- (refRatioData * totalCountData)
absData[altCountData<10] <- NA
absData[refCountData<10] <- NA

# Not necessary: remove data dat is not needed anymore
rm(list=setdiff(ls(), c("absData","clindata","annData","DATA.DIR","RES.DIR","totalCountData", "genes_of_interest")))

# Filter out SNPs with only NAs in absData (i.e. not measured or only homozygotes)
selection <- which(apply(absData, 1, function(x) all(is.na(x))))
totalCountData.filt <- totalCountData[-selection, ]
absData.filt <- absData[-selection, ]

# remove X-chromosome
absData.filtx <- absData.filt[!grepl("^X:", rownames(absData.filt)), ]
totalCountData.filtx <- totalCountData.filt[!grepl("^X:", rownames(totalCountData.filt)), ]


# Calculate the median ASE value among all samples for binomial model
absMedians <- apply(absData.filtx, 2, median, na.rm = T)
hist(absMedians)
IS.MEDIAN = median(absMedians)

# Export all SNPs for future reference
setwd(RES.DIR)
write.table(rownames(absData.filtx), file = "allAseSnps.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

# Here we need to load the position to rsId mapping from external file
# see instructions here: https://github.com/ashviyer/ASE-network-clustering/blob/main/Data/readme-map-pos2rsid.txt
# rename the rownames of absdata with rsids
setwd(DATA.DIR)
absData.filtx$positions <- rownames(absData.filtx)
mapping.rsids <- read.delim("mapping-rsids.tsv")
merged_absData <- merge(absData.filtx, mapping.rsids, by.x = "positions", by.y ="pos")
merged_absData <- merged_absData %>% filter(!is.na(rsid), rsid != "", rsid != "NA")
rownames(merged_absData) <- merged_absData$rsid
merged_absData$positions <- NULL
merged_absData$rsid <- NULL

merged_totalCountData <- merge(totalCountData.filtx, mapping.rsids, by.x = 0, by.y ="pos")
merged_totalCountData <- merged_totalCountData %>% filter(!is.na(rsid), rsid != "", rsid != "NA")
rownames(merged_totalCountData) <- merged_totalCountData$rsid
merged_totalCountData$rsid <- NULL
merged_totalCountData$Row.names <- NULL

#-----------------------------------------------------------------------------#
# 3. From SNPs to genes: functions and annotation
#-----------------------------------------------------------------------------#

# Create giant table, once, linking all ASE SNPs to all genes
#-----------------------------------------------------------------------------#

# IMPORTANT: THIS PART NEEDS TO BE RERUN FOR EVERY CHANGE IN INCLUDED SAMPLES
grch37.snp <<- useMart(biomart = "ENSEMBL_MART_SNP", 
                         dataset = "hsapiens_snp", host = "https://grch37.ensembl.org")
grch37.gene <<- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                          dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
bmSnpData <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"), 
                     filters = "snp_filter", 
                     values = rownames(merged_absData),
                     mart = grch37.snp)
  

bmGeneData <- getBM(attributes = c("hgnc_symbol", "external_gene_name", 
                                     "ensembl_gene_id","description"), 
                      filters = "ensembl_gene_id", 
                      values = unique(bmSnpData$ensembl_gene_stable_id), 
                      mart = grch37.gene)
bmGeneData[bmGeneData[, 1] == "", 1] <- bmGeneData[bmGeneData[, 1] == "", 2]
colnames(bmSnpData)[2] <- "ensembl_gene_id"
snp2geneData <- merge(bmSnpData, bmGeneData[, c(1, 3)], 
                        by = "ensembl_gene_id")

colnames(snp2geneData)[2] <- "rsid"

setwd(DATA.DIR)
save(list = "snp2geneData", file = "snp2geneData.RData")

load("snp2geneData.RData")

#-----------------------------------------------------------------------------#
# 4. popupation median based metric
#-----------------------------------------------------------------------------#
# Population median as expected!

# Binomial test for significance of ASE columns are samples, rows are rsIDs
ase.binom.test <- function(countdata, asedata, expected) {
  asePvalueMatrix <- c()
  for (i in colnames(countdata)) {
    nr  <- countdata[, i]
    absind <- asedata[, i]
    p <- apply(cbind(absind, nr), 1, function(r) {
      if (is.na(r[1])) {
        NA
      } else {
        binom.test(round(r[1] * r[2]), r[2], 
                   expected)$p.value
      }
    })
    asePvalueMatrix <- cbind(asePvalueMatrix, p)
    paste0(i)
  }
  rownames(asePvalueMatrix) <- rownames(asedata)
  colnames(asePvalueMatrix) <- colnames(asedata)
  return(asePvalueMatrix)
}

asePvalueMatrix <- ase.binom.test(countdata = merged_totalCountData,
                                  asedata = merged_absData,
                                  expected = IS.MEDIAN)

# Correct p-values
aseQvalueMatrix <- asePvalueMatrix
aseQvalueMatrix[!is.na(aseQvalueMatrix)] <-  qvalue(as.vector(as.matrix(aseQvalueMatrix[!is.na(aseQvalueMatrix)])))$qvalues

# Remove p- and q-values for loci with ASE < Median
aseQvalueMatrix <- as.matrix(aseQvalueMatrix)
asePvalueMatrix <- as.matrix(asePvalueMatrix)

merged_absData <- as.matrix(merged_absData)

aseQvalueMatrix[merged_absData<IS.MEDIAN] <- NA
asePvalueMatrix[merged_absData<IS.MEDIAN] <- NA

aseQvalueMatrix <- as.data.frame(aseQvalueMatrix)
asePvalueMatrix <- as.data.frame(asePvalueMatrix)

merged_absData <- as.data.frame(merged_absData)

# Check distributions of ASE and p-values
absScores <- merged_absData
absScores$rsid <- rownames(absScores)
absScores <- pivot_longer(absScores, c(1:(ncol(absScores)-1)), values_to = "ASE", names_to = "id")
absScores <- na.omit(absScores)
absScores <- full_join(absScores,snp2geneData)

# Density plot for ASE values in all genes
densitytotal <- ggplot(absScores, aes(x = ASE)) +
  geom_density(linewidth = 0.8, alpha = 0.2) +
  labs(y = "Density", x = "ASE") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  )
tiff("ASE_Density_Total.tif", res = 300, width = 85, height = 85, units = "mm")
print(densitytotal)
dev.off()

# Create matrix and merge with gene annotations
apm <- as.data.frame(asePvalueMatrix)
apm$rsid <- rownames(apm)
apm <- pivot_longer(apm, c(1:(ncol(apm)-1)),names_to = "id")
colnames(apm)[3] <- "pvalue"
apm <- na.omit(apm)
snps2genesP <- full_join(apm,snp2geneData)
snps2genesP <- snps2genesP[!is.na(snps2genesP$pvalue),]

snps2genesP$logp <- (-log10(snps2genesP$pvalue))

setwd(RES.DIR)

# Link all q-values to ENSGID for later analyses
aqm <- as.data.frame(aseQvalueMatrix)
aqm$rsid <- rownames(aqm)
aqm <- pivot_longer(aqm, c(1:(ncol(aqm)-1)),names_to = "id")
colnames(aqm)[3] <- "qvalue"
aqm <- na.omit(aqm)
aqm$qvalue[aqm$qvalue==0] <- min(aqm$qvalue[!aqm$qvalue==0])
snps2genesQ <- full_join(aqm,snp2geneData)

# For individual q-value ASE Manhattan plot
# Get annotation data necessary for plot
annData <- annData[,c(1,2,3)]
colnames(annData)[3] <- "pos"

annData <- merge(annData, mapping.rsids, by.x = "pos", by.y ="pos")
annData <- annData %>% filter(!is.na(rsid), rsid != "", rsid != "NA")
#annData[,1] <- gsub("chr", "", annData[,1])
#annData[,1] <- gsub("_.*", "", annData[,1])
annData <- annData[annData$rsid %in% rownames(aseQvalueMatrix),]
annData$chromosome <- as.integer(annData$chromosome)
annData$position <- as.integer(annData$position)
annData <- full_join(annData,snp2geneData)

setwd(RES.DIR)

#-----------------------------------------------------------------------------#
# 5. Translate ASE to 0 and 1
#-----------------------------------------------------------------------------#


absQValues <- aseQvalueMatrix
absQValues$rsid <- rownames(absQValues)
absQValues <- pivot_longer(absQValues, c(1:(ncol(absQValues)-1)), values_to = "qValue", names_to = "id")
absQValues <- na.omit(absQValues)
absQValues <- full_join(absQValues,snp2geneData)
absQValues <- na.omit(absQValues)

complete.ase.data <- merge(absQValues, absScores)

complete.ase.data <- complete.ase.data[order(complete.ase.data$id, complete.ase.data$ensembl_gene_id, complete.ase.data$ASE, decreasing = T),]
#complete.ase.data <- complete.ase.data[order(complete.ase.data$id, complete.ase.data$ensembl_gene_id, complete.ase.data$qValue, decreasing = F),]

# Take most significant hit per person per gene
complete.ase.data.filt <- complete.ase.data %>% distinct(id, ensembl_gene_id, .keep_all= TRUE)

complete.ase.data.filt$imbalanced[complete.ase.data.filt$qValue < 0.05] <- 1
complete.ase.data.filt$imbalanced[!complete.ase.data.filt$qValue < 0.05] <- 0

result <- dcast(complete.ase.data.filt, id ~ imbalanced, fun.aggregate = length, value.var = "imbalanced")
colnames(result) <- c("id", "count_0", "count_1")

plot_data <- result %>%
  pivot_longer(cols = c(count_0, count_1), 
               names_to = "value", 
               values_to = "count")

ggplot(plot_data, aes(x = id, y = count, fill = value)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("count_0" = "lightblue", "count_1" = "coral"),
                    labels = c("0", "1")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of 0s and 1s per Subject",
       x = "Subject ID", y = "Proportion", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# distribution of ASE scores that are considered imbalanced now
imbalanced.genes <- complete.ase.data.filt[complete.ase.data.filt$imbalanced == 1,]
hist(imbalanced.genes$ASE)

save(list = "complete.ase.data.filt", file = "complete.ase.data.filt.RData")


#load("complete.ase.data.filt.RData")





############# CHECK ALL THIS CODE BELOW - NEEDED? ##############################

# 
# absGene <- absScores[,c("id","ASE","ensembl_gene_id")]
# absGene <- na.omit(absGene)
# absGene <- absGene[order(absGene$id, 
#                          absGene$ensembl_gene_id, 
#                          absGene$ASE, decreasing = T),]#add pvalue here
# 
# # Take most significant hit per person per gene
# absGene <- absGene %>% distinct(id, ensembl_gene_id, .keep_all= TRUE)
# absGene <- na.omit(absGene)
# # Check if clustering is possible just on ASE scores (NA = 0)
# absGene <- pivot_wider(absGene, names_from = id, values_from = ASE)
# absGene <- as.data.frame(absGene)
# rownames(absGene) <- absGene$ensembl_gene_id
# absGene <- absGene[,-1]
# absGene <- t(absGene)
# absGene[is.na(absGene)] <- IS.MEDIAN
# 
# # install.packages("factoextra")
# library(factoextra)
# 
# library(readxl)
# library(tidyverse)
# library(ggfortify)
# library(ggplot2)
# library(factoextra)
# library(cluster)
# 
# ##distance matrix
# distance <- get_dist(absGene, stand = TRUE, method = "euclidean")
# 
# png("distance_matrix_ASE.png", width = 1500, height = 1000)
# fviz_dist(distance, order=TRUE, gradient = list(mid = "white", high = "#FC4E07")) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))
# dev.off()
# 
# ##No apparent clusters
# 
# ##Determining the number of optimal clusters
# set.seed(123)
# 
# fviz_nbclust(absGene,kmeans, method = "wss", k.max = 10, nboot = 100)
# fviz_nbclust(absGene,kmeans, method = "silhouette", k.max = 10, nboot = 100)
# 
# ##gap statistics for kmeans
# gap_stat <- clusGap(absGene, FUN = kmeans, nstart = 28, K.max = 10, B = 200)
# 
# png("gap_plot_ASE.png", width = 1500, height = 1000)
# fviz_gap_stat(gap_stat) +
#   theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))
# dev.off()
# 
# # Gap statistic for hierarchical clustering
# gap_stat_hclust <- clusGap(absGene, FUN = hcut, K.max = 10, B = 100)
# fviz_gap_stat(gap_stat_hclust)
# 
# Take table that keeps only max ASE for multiple SNPs in same Gene
# Take genes that pass filtering of RNA-seq quality control




setwd(DATA.DIR)
##keeping the protein codign genes from the count data
adjusted_counts<- read.delim("protein_coding_genes.txt")
kept_genes<- adjusted_counts$x




############# CHECK ALL THIS CODE BELOW - NEEDED? ##############################

#cv_genes <- read.table("cv_genes.txt")$gene
#diffgenes <- read.table("gene_list_diffase.txt", header = T)
#diffgenes <- unique(diffgenes$x)
#table(diffgenes %in% cv_genes)
#imbalance <- read.table("mutationmatrix_900_cvgenes.txt", sep = "\t")


#netgenes <- unique(imbalance$V2)
netgenes <- unique(imbalanced.genes$ensembl_gene_id)


#table(netgenes %in% cv_genes)


#kept_genes <- c(kept_genes, cv_genes)

networkdata <- full_join(absScores, snps2genesQ)
networkdata <- networkdata[networkdata$ASE>=IS.MEDIAN,]
networkdata <- networkdata[networkdata$ensembl_gene_id %in% kept_genes,]
genesQsample <- networkdata[,c("id","qvalue","ensembl_gene_id","ASE")]
genesQsample <- genesQsample[order(genesQsample$id,
                                   genesQsample$ensembl_gene_id,
                                   genesQsample$qvalue),]

# Take most significant hit per person per gene
genesQsample <- genesQsample %>% distinct(id, ensembl_gene_id, .keep_all= TRUE)
genesQsample <- na.omit(genesQsample)


setwd(DATA.DIR)

networkdata <- genesQsample

# Set significant imbalance to 1 (mutated), rest to 0 (not mutated)
networkdata$mutation <- NA
networkdata$mutation[networkdata$qvalue<0.05] <- 1
networkdata$mutation[!networkdata$qvalue<0.05] <- 0

networkdata <- networkdata[,c("id", "ensembl_gene_id", "mutation")]
networkdata <- networkdata[networkdata$mutation==1,]

networkdata <- networkdata[networkdata$ensembl_gene_id %in% kept_genes,]

# 1. getSTRINGdb for human
# Apply 0.9 confidence score
string_db <- STRINGdb$new(species=9606, version = "12", score_threshold = 900)
human_graph <- string_db$get_graph()

# 2. create adjacency matrix
adj_matrix <- igraph::as_adjacency_matrix(human_graph)


# 3. map gene ids to protein ids
### extract protein ids from the human network
protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                      function(x) x[2])

### get protein to gene id mappings
mart_results <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                      filters = "ensembl_peptide_id", 
                      values = protein_ids,
                      mart = grch37.gene)

### replace protein ids with gene ids
ix <- match(protein_ids, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]

newnames <- protein_ids
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
  mart_results[ix, 'ensembl_gene_id']
rownames(adj_matrix) <- newnames
colnames(adj_matrix) <- newnames

ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
nullrows <- Matrix::rowSums(ppi)==0
ppi <- ppi[!nullrows,!nullrows] # ppi is the network with gene ids

# CHANGE TO DIRECTED WITH TOY STAR-LIKE NETWORK
graph <- graph_from_adjacency_matrix(ppi, mode = "undirected")
edge_list <- as.data.frame(as_edgelist(graph))

library(Rgraphviz)
library(CePa)

setwd(DATA.DIR)

# Filter ppi on only genes accepted in RNAseq qc
edge_list <- edge_list[edge_list$V1 %in% kept_genes,]
edge_list <- edge_list[edge_list$V2 %in% kept_genes,]
edge_list <- unique(edge_list)
edge_list <- as.matrix(edge_list)

graph_final <- graph_from_edgelist(as.matrix(edge_list))
components <- igraph::components(graph_final)
biggest_cluster_id <- which.max(components$csize)

# ids
vert_ids <- V(graph_final)[components$membership == biggest_cluster_id]

# subgraph
graph_reduced <- induced_subgraph(graph_final, vert_ids)
edge_list_reduced <- as_edgelist(graph_reduced)

geneids <- c(edge_list_reduced[,1], edge_list_reduced[,2])
geneids <- unique(geneids)##unique genes in the componenet

networkdata <- networkdata[networkdata$ensembl_gene_id %in% geneids,]

write.table(networkdata[,c(1,2)],"mutationmatrix_900_wasp.txt", quote = F, 
            col.names = F, row.names = F, sep = "\t")
write.table(edge_list_reduced,"ppi_network_900_wasp.txt", col.names = F, row.names = F, quote = F, sep = "\t")


##plots for the genes selected in the componenet of the ppi network
#histogram of medians of selected genes
#calculate median value for ASE/sample
significant_ase <- genesQsample[genesQsample$ensembl_gene_id %in% geneids,]
  
median_ase_by_id <- aggregate(ASE ~ id, data = significant_ase, FUN = median, na.rm = TRUE)
hist(median_ase_by_id$ASE, main = "Histogram for the median significant ASE scores persample")

significant_ase$imbalanced <- ifelse(significant_ase$qvalue < 0.05, 1L, 0L)
significant_ase_only1<-subset(significant_ase, imbalanced==1)
hist(significant_ase_only1$ASE)





result <- dcast(significant_ase, id ~ imbalanced, fun.aggregate = length, value.var = "imbalanced")
colnames(result) <- c("id", "count_0", "count_1")

plot_data <- result %>%
  pivot_longer(cols = c(count_0, count_1), 
               names_to = "value", 
               values_to = "count")

ggplot(plot_data, aes(x = id, y = count, fill = value)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("count_0" = "lightblue", "count_1" = "coral"),
                    labels = c("0", "1")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of 0s and 1s per Subject",
       x = "Subject ID", y = "Proportion", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# df1$colA and df2$colB are character columns
x <- trimws(as.character(mutationmatrix_900_wasp$V2))
y <- trimws(as.character(ppi_network_900_wasp$V1))


# 2) Number of unique matching values (intersection size)
n_match_unique <- length(intersect(unique(na.omit(x)), unique(na.omit(y))))
#-----------------------------------------------------------------------------#
# 6. Identify differential imbalance between groups
#-----------------------------------------------------------------------------#

setwd(DATA.DIR)
clindata <- "HTT_900_notwasp_3_cluster_assignments.csv"
clindata <- read.table(clindata, sep = ",", header = F)
colnames(clindata) <- c("id", "cluster")
clindata$cluster <- as.factor(clindata$cluster)

absScores_merged <- full_join(absScores,clindata)
absScores_merged<- absScores_merged[complete.cases(absScores_merged),]

# Density plot for ASE values in all genes
densitytotal <- ggplot(absScores_merged, aes(x = ASE, fill = cluster)) +
  geom_density(linewidth = 0.5, alpha = 0.5) +
  labs(y = "Density", x = "ASE") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#DB5F57","#57D3DB","pink" )) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  )

tiff("ASE_Density_3.tif", res = 300, width = 85, height = 85, units = "mm")
print(densitytotal)
dev.off()

metadata <- read.table("clinical.data.final.txt", header = T, sep = "\t")
metadata <- metadata[,c("Run", "Dataset.batch")]
colnames(metadata) <- c("id","batch")
metadata <- metadata[metadata$id %in% absScores_merged$id,]

#absScores <- full_join(absScores,clindata)
#absScores$cluster <- as.factor(absScores$cluster)

absScores_merged <- full_join(absScores_merged,metadata)
absScores_merged$batch <- as.factor(absScores_merged$batch)

aqmCluster <- full_join(clindata,aqm)
aqmCluster <- aqmCluster[aqmCluster$qvalue<0.05,]
aqmCluster1 <- aqmCluster[aqmCluster$cluster==1,]
aqmCluster2 <- aqmCluster[aqmCluster$cluster==2,]
aqmCluster3 <- aqmCluster[aqmCluster$cluster==3,]

# Random checks of number of measured loci per cluster
check1 <- as.data.frame(table(aqmCluster1$id))
check2 <- as.data.frame(table(aqmCluster2$id))
check3 <- as.data.frame(table(aqmCluster3$id))

median(check1$Freq)
median(check2$Freq)
median(check3$Freq)

absScoresMed <- absScores_merged[absScores_merged$ASE>=IS.MEDIAN,]
absScoresMed <- full_join(clindata,absScoresMed)
absScores1 <- absScoresMed[absScoresMed$cluster==1,]
absScores2 <- absScoresMed[absScoresMed$cluster==2,]
absScores3 <- absScoresMed[absScoresMed$cluster==3,]

check1 <- as.data.frame(table(absScores1$id))
check2 <- as.data.frame(table(absScores2$id))
check3 <- as.data.frame(table(absScores3$id))

median(check1$Freq)
median(check2$Freq)
median(check3$Freq)

check <- as.data.frame(table(absScoresMed$id))
median (check$Freq)

##perform non parametric anova across the three clusters
# ========================================
# STEP 1: Prepare Data
# ========================================

# Since rsid already exists, don't use rownames_to_column
ase_long <- absScores_merged %>%
  filter(!is.na(ASE)) %>%  # Remove missing values
  mutate(cluster = factor(cluster, levels = c("1", "2", "3")))

# Check the result
head(ase_long)
dim(ase_long)

# Check how many clusters each SNP appears in
snp_cluster_counts <- ase_long %>%
  group_by(rsid) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    clusters_present = paste(sort(unique(as.character(cluster))), collapse = ",")
  )

filter_cluster_snp<-subset(snp_cluster_counts, n_clusters ==3)

ase_long <- absScores_merged %>%
  filter(!is.na(ASE)) %>%  # Remove missing values
  mutate(cluster = factor(cluster, levels = c("1", "2", "3")))

# Check the result
head(ase_long)
dim(ase_long)

# Check how many clusters each SNP appears in
snp_cluster_counts <- ase_long %>%
  group_by(rsid) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    clusters_present = paste(sort(unique(as.character(cluster))), collapse = ",")
  )

filter_cluster_snp<-subset(snp_cluster_counts, n_clusters ==3)


ase_long_1 <- ase_long %>%
  semi_join(filter_cluster_snp %>% distinct(rsid), by = "rsid")

# ========================================
# STEP 2: Kruskal-Wallis Test per SNP
# ========================================

cat("\nPerforming Kruskal-Wallis test for each SNP...\n")

# Perform Kruskal-Wallis test for each rsid
kw_results <- ase_long %>%
  group_by(rsid) %>%
  do(tidy(kruskal.test(ASE ~ cluster, data = .))) %>%
  ungroup() %>%
  mutate(
    padj_BH = p.adjust(p.value, method = "BH"),
    padj_bonferroni = p.adjust(p.value, method = "bonferroni")
  ) %>%
  arrange(p.value)

# View results
cat("\nTop 20 results:\n")
print(head(kw_results, 20))

# Function for differential ASE analysis BETWEEN 2 GROUPS
diffAse <- function(asedata, clinicaldata, groupVar) {
  gvar <- data.frame(clinicaldata[, c("id", groupVar)])
  gvar2 <- gvar[match(colnames(asedata), gvar[, 1]),]
  pVals <- as.data.frame(apply(asedata, 1, function(x) {
    testData <- split(as.vector(as.matrix(x)), gvar2[, groupVar])
    p <- NA
    try ({
      p <- wilcox.test(testData[[1]], 
                       testData[[2]], na.rm = T)$p.value
    }, silent = TRUE
    )
    for (j in unique(gvar[, groupVar])) {
      if (length(na.omit(testData[[j]]))<3) {
        p <- NA
      } else {
        p <- p
      }
    }
    return(p)
  }))
  colnames(pVals) <- "pvalue"
  pVals$rsid <- rownames(pVals)
  return(pVals)
}

# Function for differential ASE effect size BETWEEN 2 GROUPS (effect is d(medians))
diffAseEffect <- function(asedata, clinicaldata, groupVar) {
  gvar <- data.frame(clinicaldata[, c("id", groupVar)])
  gvar2 <- gvar[match(colnames(asedata), gvar[, 1]),]
  eVals <- as.data.frame(apply(asedata, 1, function(x) {
    testData <- split(as.vector(as.matrix(x)), gvar2[, groupVar])
    e <- NA
    try ({
      e <- (median(testData[[2]], na.rm = T) - median(testData[[1]], 
                                                      na.rm = T))
    }, silent = TRUE
    )
    for (j in unique(gvar[, groupVar])) {
      if (length(na.omit(testData[[j]]))<3) {
        e <- NA
      } else {
        e <- e
      }
    }
    return(e)
  }))
  colnames(eVals) <- "effect"
  eVals$rsid <- rownames(asedata)
  return(eVals)
}

diffAseRes <- diffAse(merged_absData, clindata, "cluster")
diffAseEff <- diffAseEffect(merged_absData, clindata, "cluster")

# Make manhattan plot for current comparison
diffAseRes <- na.omit(diffAseRes)
annTest <- annData[annData$rsid %in% diffAseRes$rsid,]

#############################
#FOR 3 CLUSTERS anova
#############################

library(dplyr)
library(tidyr)
library(broom)

kw_results <- ase_long %>%   # << use the filtered data
  group_by(rsid) %>%
  do(broom::tidy(kruskal.test(ASE ~ cluster, data = .))) %>%
  ungroup() %>%
  mutate(
    padj_BH = p.adjust(p.value, method = "BH"),
    padj_bonferroni = p.adjust(p.value, method = "bonferroni")
  ) %>%
  arrange(padj_BH, p.value)

eps2_kw <- function(H, n, k) {
  e <- (as.numeric(H) - k + 1) / (n - k)
  ifelse(is.finite(e) & e > 0, e, 0)
}

diffAseKW3 <- function(asedata, clinicaldata, groupVar = "cluster",
                       min_n_per_group = 3,
                       levels = c("1","2","3")) {
  
  asedata <- as.matrix(asedata)
  
  gvar <- clinicaldata[, c("id", groupVar)]
  idx <- match(colnames(asedata), gvar$id)
  
  if (all(is.na(idx))) {
    stop("No clinicaldata$id matches colnames(asedata). Check sample IDs.")
  }
  
  gvar2 <- gvar[idx, , drop = FALSE]
  
  keep <- !is.na(gvar2[[groupVar]])
  asedata2 <- asedata[, keep, drop = FALSE]
  grp <- factor(gvar2[[groupVar]][keep], levels = levels)
  
  out <- t(apply(asedata2, 1, function(x) {
    x <- as.numeric(x)
    ok <- !is.na(x) & !is.na(grp)
    x <- x[ok]
    g <- droplevels(grp[ok])
    
    tab <- table(g)
    
    # require all 3 clusters + enough samples each
    if (length(tab) < 3 || any(tab < min_n_per_group)) {
      n_c1 <- unname(tab["1"]); if (is.na(n_c1)) n_c1 <- 0L
      n_c2 <- unname(tab["2"]); if (is.na(n_c2)) n_c2 <- 0L
      n_c3 <- unname(tab["3"]); if (is.na(n_c3)) n_c3 <- 0L
      return(c(H = NA_real_, pvalue = NA_real_, epsilon2 = NA_real_,
               n_total = length(x),
               n_c1 = as.integer(n_c1), n_c2 = as.integer(n_c2), n_c3 = as.integer(n_c3),
               med_c1 = NA_real_, med_c2 = NA_real_, med_c3 = NA_real_))
    }
    
    kt <- kruskal.test(x ~ g)
    H <- unname(kt$statistic)
    p <- kt$p.value
    n <- length(x)
    k <- nlevels(g)
    e2 <- eps2_kw(H, n, k)
    
    meds <- tapply(x, g, median, na.rm = TRUE)
    med_c1 <- unname(meds["1"])
    med_c2 <- unname(meds["2"])
    med_c3 <- unname(meds["3"])
    
    n_c1 <- as.integer(unname(tab["1"]))
    n_c2 <- as.integer(unname(tab["2"]))
    n_c3 <- as.integer(unname(tab["3"]))
    
    c(H = H, pvalue = p, epsilon2 = e2,
      n_total = n,
      n_c1 = n_c1, n_c2 = n_c2, n_c3 = n_c3,
      med_c1 = as.numeric(med_c1), med_c2 = as.numeric(med_c2), med_c3 = as.numeric(med_c3))
  }))
  
  out <- as.data.frame(out)
  out$rsid <- rownames(out)
  
  out$padj_BH <- p.adjust(out$pvalue, method = "BH")
  out$padj_bonferroni <- p.adjust(out$pvalue, method = "bonferroni")
  
  out[order(out$padj_BH, out$pvalue), ]
}

kw3_res <- diffAseKW3(merged_absData, clindata, groupVar = "cluster",
                      min_n_per_group = 3, levels = c("1","2","3"))

head(kw3_res, 20)

sig_snps <- subset(kw3_res, pvalue < 0.05)

diffAseRes <- na.omit(kw3_res)
annTest <- annData[annData$rsid %in% diffAseRes$rsid,]

#######################################################################





# Function to make Manhattan plot
# Expects data frame of rsid SID and pvalue/qvalue
# Expects data frame of chromosome, position, rsid, hgnc_symbol
# Standard statvalue = "pvalue" can be changed to "qvalue"
aseManhattan <- function(input, 
                         statvalue = "pvalue", 
                         annotation, 
                         gws, 
                         filename) {
  if (statvalue == "qvalue") {
    input <- full_join(input, annotation, by = "rsid")
    input <- input %>% arrange(chromosome, position)
    
    # Create cumulative positions for ordering x-axis
    input$BPcum <- NA
    input$BPcum <- as.numeric(input$BPcum)
    s <- 0
    nbp <- c()
    chrnum <- unique(input$chromosome)
    for (i in 1:length(chrnum)) {
      nbp[i] <- max(input[input$chromosome == chrnum[i],]$position)
      input[input$chromosome == chrnum[i],"BPcum"] <- input[input$chromosome == chrnum[i],"position"] + s
      s <- s + nbp[i]
    }
    
    #Create centering for each chromosome x-axis location
    axis.set <- input %>%
      group_by(chromosome) %>%
      summarize(center = (max(BPcum) + min(BPcum)) / 2)
    ylim <- (-log10(min(input$qvalue)) + 1)
    sig = gws
    
    # ggplot
    manhattanplot <- ggplot(input, aes(x = BPcum, y = -log10(P), 
                                       colour = CHR)) +
      geom_point(size = 1) +
      geom_point(data = input, aes(x = BPcum, y = -log10(P)), size = 1) +
      scale_x_continuous(expand = c(0.05,0.05), 
                         breaks = axis.set$center, labels = axis.set$CHR,
                         guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_color_manual(values = rep(c("#276EBF", "#183059"), 24)) + 
      geom_point(data = input[input$P<1e-5,], color="#ff9933", size = 1) +
      labs(x = "Chromosome", y = expression(paste("-log"[10],"(", plain(p),")")), title = paste0(unique(input$Analysis))) +
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.text.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.ticks = element_blank()
      )
    tiff(filename = paste0(filename,".tif"), res = 300, width = 175, height = 85, units = "mm")
    print(manhattanplot)  
    dev.off()
    
  } else {
    input <- full_join(input, annotation, by = "rsid")
    input <- input %>% arrange(chromosome, position)
    
    # Create cumulative positions for ordering x-axis
    input$BPcum <- NA
    input$BPcum <- as.numeric(input$BPcum)
    s <- 0
    nbp <- c()
    chrnum <- unique(input$chromosome)
    for (i in 1:length(chrnum)) {
      nbp[i] <- max(input[input$chromosome == chrnum[i],]$position)
      input[input$chromosome == chrnum[i],"BPcum"] <- input[input$chromosome == chrnum[i],"position"] + s
      s <- s + nbp[i]
    }
    
    #Create centering for each chromosome x-axis location
    axis.set <- input %>%
      group_by(chromosome) %>%
      summarize(center = (max(BPcum) + min(BPcum)) / 2)
    ylim <- (-log10(min(input$pvalue)) + 1)
    sig = gws
    
    # ggplot
    manhattanplot <- ggplot(input, aes(x = BPcum, y = -log10(P), 
                                       colour = CHR)) +
      geom_point(size = 1) +
      geom_point(data = input, aes(x = BPcum, y = -log10(P)), size = 1) +
      scale_x_continuous(expand = c(0.05,0.05), 
                         breaks = axis.set$center, labels = axis.set$CHR,
                         guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
      scale_color_manual(values = rep(c("#276EBF", "#183059"), 24)) + 
      geom_point(data = input[input$P<1e-5,], color="#ff9933", size = 1) +
      labs(x = "Chromosome", y = expression(paste("-log"[10],"(", plain(p),")")), title = paste0(unique(input$Analysis))) +
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.text.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.ticks = element_blank()
      )
    tiff(filename = paste0(filename,".tif"), res = 300, width = 175, height = 85, units = "mm")
    print(manhattanplot)  
    dev.off()
  }
}

diffAseRes$P <- diffAseRes$pvalue
annTest$CHR<-annTest$chromosome

# Get cummulative position
aseManhattan(input = diffAseRes,
             statvalue = "pvalue",
             annotation = annTest,
             gws = 0.05,
             filename = paste0("1_vs_2_Manhattan"))

diffAseRes <- full_join(diffAseRes, annTest)
diffAseRes <- na.omit(diffAseRes)

write.csv(diffAseRes, "diffASE_results.csv", quote = F, sep = "/t")

# Make boxplot for most significant SNP hit in across Clusters comparison
# Get necessary data for boxplot
boxplotData <- merged_absData
boxplotData$rsid <- rownames(boxplotData)
boxplotData <- pivot_longer(boxplotData, c(1:(ncol(boxplotData)-1)),
                            names_to = "id")
colnames(boxplotData)[3] <- "ASE"
colnames(boxplotData)[2] <- "id"
boxplotData <- full_join(boxplotData, clindata[,c("id","cluster")])
boxplotData <- na.omit(boxplotData)

# Get most significant hit
snpBoxplot <- 
  diffAseRes[diffAseRes$pvalue==min(diffAseRes$pvalue),"rsid"]
snpBoxplotData <- boxplotData[boxplotData$rsid==snpBoxplot,]
snpBoxplotData$cluster <- as.factor(snpBoxplotData$cluster)

summ <- snpBoxplotData %>% group_by(cluster) %>% dplyr::summarize(n = n(), avg = max(ASE))

# Create ggplot
boxplotplot <- ggplot(snpBoxplotData, aes(x = cluster, y = ASE, fill = cluster)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Cluster", y = "ASE score", fill = "Cluster") +
  scale_fill_brewer(palette = "Dark2") +
  coord_fixed(ratio = 8) +
  geom_text(data = summ, inherit.aes = FALSE, aes(x = cluster, label = n, y = avg+0.05), size = 3) +
  theme(
    legend.position = "right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  ) +
  geom_jitter(shape = 16, alpha = 0.5, width = 0.25)
tiff(filename=paste0("Boxplot_", snpBoxplot, ".tif"), 
     width=85, height=85, units = "mm", res = 300)
print(boxplotplot)
dev.off()



#-----------------------------------------------------------------------------#
# 7. Look for genes measured in at least 50% samples 
#     with a significant imbalance in at least 80%
#-----------------------------------------------------------------------------#
genesQsample <- snps2genesQ##rename the genesQsample as it already defined as something different above

genesQsample <- genesQsample[order(genesQsample$rsid),]

genesQsample <- genesQsample[genesQsample$qvalue<0.05,]
genesQsample <- genesQsample[!is.na(genesQsample$ensembl_gene_id),]

# number of unique samples
N <- length(unique(genesQsample$id))   # use your sample column name (id)

# 50% threshold
min_samples <- ceiling(0.5 * N)

snpfrq <- as.data.frame(table(genesQsample$rsid))
colnames(snpfrq) <- c("rsid","freq")

genesQsample <- dplyr::left_join(genesQsample, snpfrq, by = "rsid")
genesQsample <- genesQsample[genesQsample$freq >= min_samples, ]

genesQsample <- genesQsample[,c("id","qvalue","ensembl_gene_id")]

genesQsample <- genesQsample[order(genesQsample$id, 
                                   genesQsample$ensembl_gene_id, 
                                   genesQsample$qvalue),]

# Take most significant hit per person per gene
genesQsample <- genesQsample %>% distinct(id, ensembl_gene_id, .keep_all= TRUE)
genesQsample <- na.omit(genesQsample)

N <- length(unique(genesQsample$id))      # number of samples
min_samples_80 <- ceiling(0.8 * N)

AseGeneCounts <- as.data.frame(table(genesQsample[genesQsample$qvalue < 0.05, "ensembl_gene_id"]))
colnames(AseGeneCounts) <- c("ensembl_gene_id", "freq")

AseGeneCounts80 <- AseGeneCounts[AseGeneCounts$freq >= min_samples_80, "ensembl_gene_id"]

# Volcano plot with effect size
volDat <- absScores[,c(1:3,6)]
volDat <- full_join(volDat, apm)
volDat <- full_join(volDat, aqm)
volDat <- na.omit(volDat)
ylim <- (-log10(min(volDat$pvalue)) + 0.5)

volDat$threshold <- factor(ifelse(volDat$pvalue<0.05,
                                   1,0))

volDat1 <- volDat[volDat$cluster==1,]
volDat2 <- volDat[volDat$cluster==3,]

# ggplot volcanoplot
volcanoplot <- ggplot(volDat2, aes(x = ASE, y = -log10(qvalue))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
  geom_point(alpha = 0.5) +
  labs(x = "ASE Score", 
       y = expression(paste("-log"[10],"(", plain(q),")"))) +
  theme_minimal() + 
  theme(legend.position = "none") +
  theme( 
    legend.position = "right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0), face = "bold"),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0), face = "bold"),
    axis.ticks = element_blank()
  )
tiff(filename="volcano_cluster2.tif", width=200, height=150, units = "mm", res = 300)
print(volcanoplot)
dev.off()
