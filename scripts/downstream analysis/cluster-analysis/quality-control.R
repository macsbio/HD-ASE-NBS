#!/usr/bin/env Rscript

##import packages
library(readxl)
library(tidyverse)
library(DESeq2)
library(ggfortify)
library(ggplot2)
library(pheatmap)
library(sva)
library(biomaRt)
library(DescTools)
library(ComplexHeatmap)
library(circlize)

# set working dir to the folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# paths relative to the script directory
data_dir <- file.path(getwd(), "Data")
out_dir  <- file.path(getwd(), "output")

# create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# (optional) stop if Data folder is missing
if (!dir.exists(data_dir)) stop("Missing Data folder: ", data_dir)


##import processed count data 
qc_count_data <- read.csv(file.path(out_dir, "qc-count-data.csv"), header = TRUE, row.names = 1)

#import clinical data
clinical_data <- read.delim(file.path(out_dir, "cleaned-metadata.txt"), header = TRUE)

##import clusters
cluster_info<- read.csv(file.path(data_dir, "HTT_900_notwasp_3_cluster_assignments.csv"), header = F)

##Filtering genes####

## 1.remove lowly expressed genes 
removed_genes<-qc_count_data[!rowMeans(qc_count_data) > 10,]## select genes with row means across samples less than 10 (filtered out)
qc_count_data.filt <- qc_count_data[rowMeans(qc_count_data) > 10,]##select genes with row means across samples greater than 10


## 2. keep protein coding genes

#select ensembl databse - H.sapiens
grch37.gene <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
filters <- listFilters(grch37.gene)#list of filters
attributes <- listAttributes(grch37.gene)#list fo attributes


##gene categories 
genes<- getBM(mart=grch37.gene, values= rownames(qc_count_data.filt),
                   filters= "ensembl_gene_id",
                   attributes= c("ensembl_gene_id", "hgnc_symbol",
                                 "gene_biotype"))

##view categories table
categories_table<-table(genes$gene_biotype)
view(categories_table)#16,528 protein coding genes

##select protein coding genes list
protein.coding.genes<-subset(genes, gene_biotype =="protein_coding")

##transformed counts for each type of gene biotype
protein.coding.genes.count<-qc_count_data.filt[protein.coding.genes$ensembl_gene_id,]

protein.coding.genes.count<- protein.coding.genes.count[,c(clinical_data$Run)]

write.table(rownames(protein.coding.genes.count),
            file = file.path(out_dir, "protein_coding_genes.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")
write.table(rownames(removed_genes),
            file = file.path(out_dir, "low_expressed_genes_removed.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")




##Sample clustering#######
#distance matrix

dist_matrix <- dist(as.matrix(t(protein.coding.genes.count)), method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "ward.D2")
#plot dendogram
plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")


##DESEQ analysis 
protein.coding.genes.count <- as.matrix(protein.coding.genes.count)
storage.mode(protein.coding.genes.count) <- 'integer'

metadata<-data.frame(sample_ids = clinical_data$Run,
                     groups = clinical_data$cluster,
                     batch= clinical_data$Dataset.batch)
rownames(metadata)<-metadata$sample_ids

metadata$groups<-as.factor(metadata$groups)
metadata$batch<-as.factor(metadata$batch)

metadata$sample_ids<-NULL

##arrange coldata and count data in same rownamess
## 1. Arrange columns
protein.coding.genes.count <- protein.coding.genes.count[, rownames(metadata)]

# filters:
# 1) total counts across all samples > 10
keep_rowsum <- rowSums(protein.coding.genes.count, na.rm = TRUE) > 10

# 2) expressed (>0) in >20% of samples
keep_prop <- rowMeans(protein.coding.genes.count > 0, na.rm = TRUE) > 0.2

# apply
protein.coding.genes.count.filt <- protein.coding.genes.count[keep_rowsum & keep_prop, , drop = FALSE]

write.table(rownames(protein.coding.genes.count.filt),
            file = file.path(out_dir, "protein_coding_genes.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")

## 3. Batch correction on filtered counts
adjusted_counts <- ComBat_seq(as.matrix(protein.coding.genes.count.filt),
                              batch = metadata$batch)

## 4. DESeq2
de_data <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                                  colData   = metadata,
                                  design    = ~ groups)

de_data <- DESeq(de_data)

## 5. Normalised + log2
normalized_counts  <- counts(de_data, normalized = TRUE)
transformed_counts <- log2(normalized_counts + 1)


###########
#QC plots
###########
#Sum of counts
png("sum of counts.png", width= 1500, height = 1000)

colors= c(rep(c("tomato"), each = c(8)),rep(c("turquoise2"), each = c(3)),rep(c("purple"), each = c(9)))
boxplot(transformed_counts ,col=colors, main="counts distribution per sample",
        xlab="Samples",ylab="vst counts",cex=0.8, outline = FALSE, las=2)
colors_names<-c(rep(c("tomato"), each = c(8)),rep(c("turquoise2"), each = c(3)),rep(c("purple"), each = c(9)))
legend("topright",c("cluster1", "cluster2", "cluster3"),fill=colors_names)
dev.off()

##boxplot colsums per cluster
colsums.cluster<-as.data.frame(colSums(transformed_counts))
colsums.cluster$Sample_ids<-rownames(colsums.cluster)
sample_match<-match(colsums.cluster$Sample_ids, cluster_info$V1)
colsums.cluster$group<-cluster_info[sample_match,]
colsums.cluster$Sample_ids<-NULL

colsums.cluster$group$V1<-NULL
colsums.cluster$group$V2<-as.factor(colsums.cluster$group$V2)


######Counts sum per cluster 
png("colsumsper cluster.png", width=700, height=700)
p<-ggplot(colsums.cluster, aes(x=group$V2, y=`colSums(transformed_counts)`, color=group$V2)) +
  geom_boxplot() +
  scale_color_manual(values=c("tomato","turquoise2","purple"))
p
dev.off()

##########################################

#Filter highly varaiable genes

##Calculate coefficient of variance 

coeff_var<-data.frame(apply(transformed_counts,1, function(x) sd(x) / mean(x) * 100))
colnames(coeff_var)<-c("cv")
hist(coeff_var$cv)##plot


write.table(transformed_counts,
            file = file.path(out_dir, "norm-transformed-count.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")

write.table(adjusted_counts,
            file = file.path(out_dir, "adjusted_raw_counts.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")

write.table(rownames(transformed_counts),
            file = file.path(out_dir, "kept-genes.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")


##Count density plot
plot(density(transformed_counts[, 1]), col = 1, ylim=c(0,1), xlim=c(-1,10),
     main="Density plot for expression values per condition 
    (Normalisation -log)",
     xlab="log2(Gene Expression values)",
     ylab="Density")
colors_density<-c(rep(c(1), each = c(8)),rep(c(2), each = c(3)),rep(c(3), each = c(9)))
for (i in 2:ncol(transformed_counts)) {
  lines(density(transformed_counts[, i]), col = colors_density[i] )
}
colors_groups<-c(1,2,3)


## 1. Prepare PCA matrix (samples x features)
pca_mat <- t(transformed_counts)                          # rows = samples, cols = features
pca_mat <- pca_mat[rowSums(pca_mat) > 0, , drop = FALSE]  # drop samples with all zeros

## 2. Run PCA
pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

## 3. Build metadata (cluster, batch, sex)

meta <- clinical_data %>%
  dplyr::select(Run, cluster, Dataset.batch) %>%
  distinct(Run, .keep_all = TRUE)

rownames(meta) <- meta$Run

## Keep only samples present in PCA matrix, in correct order
meta <- meta[rownames(pca_mat), , drop = FALSE]

## Convert to factors
meta$cluster <- as.factor(meta$cluster)
meta$Dataset.batch  <- as.factor(meta$Dataset.batch)

## 4. Helper function to plot PCA for a given variable and save to PNG
plot_pca_by <- function(pca_res, meta, var, filename, colors = NULL) {
  var <- rlang::ensym(var)
  
  p <- autoplot(
    pca_res,
    data   = meta,
    colour = rlang::as_string(var),
    size   = 2,
    alpha  = 0.8,
    frame  = TRUE
  ) +
    labs(
      title = paste("PCA colored by", rlang::as_string(var)),
      x = "PC1",
      y = "PC2",
      colour = rlang::as_string(var)
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  png(filename, width = 1500, height = 1000, res = 150)
  print(p)
  dev.off()
}

## 5. Make the three separate plots

## (a) PCA by cluster
## If you really only have two clusters:
cluster_cols <- c("tomato", "blue", "purple")
plot_pca_by(
  pca_res,
  meta,
  var      = cluster,
  filename = "pca_cluster.png",
  colors   = cluster_cols
)

## (b) PCA by batch (let ggplot choose colors or define your own)
plot_pca_by(
  pca_res,
  meta,
  var      = Dataset.batch,
  filename = "pca_batch.png"
)


###############
#Heatmaps for gene expression data
######################

##plot heatmap
cluster_pca<- subset(clinical_data, select= c(Run, cluster))
rownames(cluster_pca)<-cluster_pca$Run
cluster_pca<-as.data.frame(cluster_pca[order(cluster_pca$cluster, decreasing = FALSE),])
rownames(cluster_pca)<-cluster_pca$Run
cluster_pca$Run<-NULL

transformed.heatmap<-as.data.frame(t(transformed_counts))

transformed.heatmap<-transformed.heatmap[rownames(cluster_pca),]


transformed.heatmap<-as.matrix(t(transformed.heatmap))
cluster_pca$cluster<-as.factor(cluster_pca$cluster)


##plot for all genes 
annoCOL <- list(cluster=c("1"="tomato", "2"="turquoise2", "3"="purple"))
#rg <- max(raw_count_data.heatmap);
png("heatmap-all.png", width=1000, height = 1000)
pheatmap(transformed.heatmap, 
         cluster_cols = FALSE, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "red"))(100),
         annotation_col = cluster_pca,
         annotation_colors = annoCOL,
         show_rownames = FALSE,
         clustering_distance_rows = "correlation"
)
dev.off()

##################
#Complex heatmap
################
cluster1.list<-subset(cluster_pca, cluster ==1)
cluster2.list<-subset(cluster_pca, cluster ==2)
cluster3.list<-subset(cluster_pca, cluster ==3)
cluster1 = transformed_counts[,rownames(cluster1.list)]
cluster2 = transformed_counts[,rownames(cluster2.list)]
cluster3 = transformed_counts[,rownames(cluster3.list)]

##color for heatmap
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
col_fun(seq(-3, 3))


cluster1_scaled<-t(apply(cluster1, 1, scale))
colnames(cluster1_scaled)<-colnames(cluster1)
cluster2_scaled<-t(apply(cluster2, 1, scale))
colnames(cluster2_scaled)<-colnames(cluster2)
cluster3_scaled<-t(apply(cluster3, 1, scale))
colnames(cluster3_scaled)<-colnames(cluster3)

ht1 = Heatmap(cluster1_scaled,name="Cluster 1", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              show_row_names = FALSE, show_row_dend = FALSE, use_raster = FALSE
              , clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")
ht2 = Heatmap(cluster2_scaled, name  = "Cluster 2", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              show_row_names = FALSE , show_row_dend = FALSE, use_raster = FALSE
              , clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")
ht3 = Heatmap(cluster3_scaled, name  = "Cluster 3", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              show_row_names = FALSE , show_row_dend = FALSE, use_raster = FALSE
              , clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")


ht1 + ht2 + ht3 




##Import data
##data for genes having an allelic imbalance
across_clusters<- read.delim("across_sig.txt", sep ="\t")

##subset significant p<0.05
across_clusters<-subset(across_clusters, pvalue <0.05)

across_clusters<-across_clusters[complete.cases(across_clusters),]

match_names <- match(across_clusters$ensembl_gene_id, rownames(transformed_counts))

ase_data <- transformed_counts[match_names,]


##plot for genes with a significant ase score
annoCOL <- list(cluster=c("1"="tomato", "2"="turquoise2"))

png("heatmap-acrossclusters.png", width=800, height=1000)
pheatmap(t(ase_data), 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("navy", "white", "red"))(100),
         annotation_col = cluster_pca,
         annotation_colors = annoCOL,
         show_rownames = FALSE,
         scale = "row",
         clustering_distance_rows = "euclidean"
)
dev.off()

cluster1.list<-subset(cluster_pca, cluster ==1)
cluster2.list<-subset(cluster_pca, cluster ==2)
cluster1 = ase_data[,rownames(cluster1.list)]
cluster2 = ase_data[,rownames(cluster2.list)]


##color for heatmap
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
col_fun(seq(-3, 3))


cluster1_scaled<-t(apply(cluster1, 1, scale))
colnames(cluster1_scaled)<-colnames(cluster1)
cluster2_scaled<-t(apply(cluster2, 1, scale))
colnames(cluster2_scaled)<-colnames(cluster2)

ht1 = Heatmap(cluster1_scaled,name="Cluster 1", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              show_row_names = FALSE, show_row_dend = FALSE, use_raster = FALSE
              , clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")
ht2 = Heatmap(cluster2_scaled, name  = "Cluster 2", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              show_row_names = FALSE , show_row_dend = FALSE, use_raster = FALSE
              , clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean")


ht1 + ht2 
