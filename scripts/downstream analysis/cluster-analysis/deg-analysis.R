
##import packages
library(readxl)

library(tidyverse)
library(ggfortify)
library(ggplot2)
library(factoextra)

library(cluster)
library(DESeq2)
#install.packages("ashr")
library(ashr)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("gprofiler2")
library(gprofiler2)
library(GOSemSim)
library(enrichplot)


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
raw.counts <- as.matrix(read.delim(file.path(out_dir, "adjusted_raw_counts.txt"), sep ="\t"))
clinical_data <- read.delim(file.path(out_dir, "cleaned-metadata.txt"), sep ="\t")



##DESEQ analyis
storage.mode(raw.counts) <- 'integer'
#create a metadata dataframne with the sample ids and other variables you would like to include in the design
metadata<-data.frame(sample_ids = clinical_data$Run,
                     groups = clinical_data$cluster,
                     batch = clinical_data$Dataset.batch)#make sure to add here 'sex' variable
rownames(metadata)<-metadata$sample_ids #add rownames to your metadata
metadata$groups<-as.factor(metadata$groups) # convert to factor
metadata$batch<-as.factor(metadata$batch) # convert to factor

metadata$sample_ids<-NULL

##arrange coldata and count data in same rownames
raw.counts<-raw.counts[,rownames(metadata)]


#deseq data object (in the design you add 'sex' as a covariate as shown below)
de_data<-DESeqDataSetFromMatrix(countData = raw.counts,
                                colData = metadata,
                                design = ~groups  
)


# Global “ANOVA” for cluster (tests whether cluster improves fit beyond batch)
dds_lrt <- DESeq(de_data, test = "LRT", reduced = ~ 1)
dds_wald <- DESeq(de_data, test = "Wald")
res_lrt <- results(dds_lrt)
res_lrt<- as.data.frame(res_lrt)
res_wald<-results(dds_wald)# pvalue is for the global 3-group cluster term
res_lrt <- res_lrt[order(res_lrt$padj), ]

sig_wald<- subset(as.data.frame(res_wald), pvalue < 0.05)

sig_genes_global <- subset(as.data.frame(res_lrt), pvalue < 0.05)

sig_genes_global$gene.id<-rownames(sig_genes_global)




#deg analysis
de_data <- DESeq(de_data)


##comparison to level 1
resultsNames(de_data)

#get the deseq output for the group comparison 
res_1vs2<-results(de_data, name="groups_2_vs_1")


##data wrangling for all comparisons
logfoldchange_2vs1<-data.frame(res_1vs2@listData)

rownames(logfoldchange_2vs1)<-res_1vs2@rownames


colnames(logfoldchange_2vs1)<-paste(colnames(logfoldchange_2vs1), "2vs1", sep = "_")


#get them in one dataframe
logfoldchange.comparison<-logfoldchange_2vs1
logfoldchange.comparison$gene.ids<-rownames(logfoldchange.comparison)

##significant deg genes
sig.diff.genes<- subset(logfoldchange.comparison, pvalue_2vs1 < 0.05)


dds_wald <- DESeq(de_data, test = "Wald")

res_2v1 <- results(dds_wald, contrast = c("groups","2","1"))
res_3v1 <- results(dds_wald, contrast = c("groups","3","1"))
res_3v2 <- results(dds_wald, contrast = c("groups","3","2"))

# optional: shrink LFC for stability
library(apeglm)
res_2v1_shr <- lfcShrink(dds_wald, contrast = c("groups","2","1"), type="apeglm")




##plot the volcano plot for all comparisons
png("volcanoplot2vs1.png", width = 1500, height = 1000)
EnhancedVolcano(logfoldchange.comparison, 
                lab = NA,
                x = 'log2FoldChange_2vs1', y = 'pvalue_2vs1', 
                pCutoff = 0.05, 
                FCcutoff = 0.5,
                ylab = bquote(~-log[10] ~ italic(p-value)),
                xlab = bquote(log[2]~ foldchange),
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                      ~ log[2] ~ FC)),
                col = c('grey30', 'grey30', 'grey30', 'red2'),
                title='Volcano plot cluster 1 vs 2',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                labSize = 7,
                ylim = c(0, 3),
                xlim = c( -3,3),
                axisLabSize = 40,
                captionLabSize = 40,
                legendLabSize = 40)
dev.off()



##import ase data 
diff.ase <- read.csv(file.path(data_dir, "diffASE_results.csv"), header = TRUE, row.names = 1)

sig.diff.ase<- subset(diff.ase, pvalue <= 0.05)


##sig expressed genes with sig allele imbalance
genes_intersect<- intersect(sig_genes_global$gene.id, sig.diff.ase$ensembl_gene_id)

clusterprof_ase<-enrichGO(gene= unique(sig.diff.ase$ensembl_gene_id),
                          keyType = "ENSEMBL",
                          ont = "BP",
                          pAdjustMethod = "fdr",
                          OrgDb = "org.Hs.eg.db",
                          universe = unique(diff.ase$ensembl_gene_id)
)

d <- godata('org.Hs.eg.db', ont="BP")
clusterprof_ase_pairwise <- pairwise_termsim(clusterprof_ase, semData = d)

png("ase-go.png", width = 1500, height = 800)
treeplot(clusterprof_ase_pairwise, showCategory = 49)
dev.off()

clusterprof_gene<-enrichGO(gene= unique(sig_genes_global$gene.id),
                          keyType = "ENSEMBL",
                          ont = "BP",
                          pAdjustMethod = "fdr",
                          OrgDb = "org.Hs.eg.db",
                          universe = unique(rownames(res_lrt))
)

clusterprof_gene_pairwise <- pairwise_termsim(clusterprof_gene, semData = d)

png("genediff-go.png", width = 1500, height = 1000)
treeplot(clusterprof_gene_pairwise, showCategory = 51)
dev.off()

clusterprof_aseandgene<-enrichGO(gene= genes_intersect,
                           keyType = "ENSEMBL",
                           ont = "BP",
                           pAdjustMethod = "none",
                           OrgDb = "org.Hs.eg.db",
                           pvalueCutoff = 0.01
)

png("aseandgene-go.png", width = 500, height = 500)
dotplot(clusterprof_aseandgene, showCategory = 49)
dev.off()

write.table(genes_intersect,
            file = file.path(out_dir, "sig.genes.ase.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t", col.names  = F)




transformed_counts<-read.delim(file.path(out_dir, "norm-transformed-count.txt"), sep ="\t")
transformed_counts<-as.data.frame(t(transformed_counts))

clusters<-subset(clinical_data, select = c(Run, cluster))
rownames(clusters)<-clusters$Run
clusters$Run<-NULL

transformed_counts$cluster<-clusters[rownames(transformed_counts),] 

#select for gens mbp and cacb
select_counts<-subset(transformed_counts, select= c("ENSG00000081026","ENSG00000234745","ENSG00000196636","ENSG00000112576","ENSG00000008018"))
select_counts$cluster<-as.factor(clusters[rownames(select_counts),])


p <- ggplot(select_counts, aes(x = cluster, y= ENSG00000081026, fill = cluster)) + 
  geom_boxplot()+
  theme_classic()+
  labs(x="Cluster", y = "MAGI3 gene expression")

png("MAGI3 expression.png", res=200)
p
dev.off()


q<-ggplot(select_counts, aes(x=cluster, y=ENSG00000234745, fill=cluster)) + 
  geom_boxplot()+
  theme_classic()+
  labs(x="Cluster", y = "HLA-B gene expression") 

png("HLA-B expression.png", res=200)
q
dev.off()


S<-ggplot(select_counts, aes(x=cluster, y=ENSG00000196636, fill=cluster)) + 
  geom_boxplot()+
  theme_classic()+
  labs(x="Cluster", y = "SDHAF3 gene expression") 

png("SDHAF3 expression.png", res=200)
S
dev.off()

t<-ggplot(select_counts, aes(x=cluster, y=ENSG00000112576, fill=cluster)) + 
  geom_boxplot()+
  theme_classic()+
  labs(x="Cluster", y = "CCND3 gene expression") 

png("CCND3 expression.png", res=200)
t
dev.off()


U<-ggplot(select_counts, aes(x=cluster, y=ENSG00000008018, fill=cluster)) + 
  geom_boxplot()+
  theme_classic()+
  labs(x="Cluster", y = "PSMB1 gene expression") 

png("PSMB1 expression.png", res=200)
U
dev.off()
