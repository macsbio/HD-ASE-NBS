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

##import data and data wrangling
raw_count_data <- read.delim(file.path(data_dir, "WASP_featureCounts-notWASP.txt"),
                             , comment.char="#",
                             header = T, row.names = 1)

raw_count_data.filt <- raw_count_data[, !names(raw_count_data) %in% c("Chr", "Start", "Strand", "End", "Length")] #remove columns 
colnames(raw_count_data.filt) <- sub(".*(SRR[0-9]+)_RG\\.bam$", "\\1", colnames(raw_count_data.filt))


#import metadata
huntington_metadata <- readxl::read_excel(file.path(data_dir, "metadata.xlsx"))

sample.match  <- read.delim(file.path(data_dir, "sample_match.txt"), sep = "\t", stringsAsFactors = FALSE)

cluster_info  <- read.csv(file.path(data_dir, "HTT_900_notwasp_3_cluster_assignments.csv"),
                          header = FALSE, stringsAsFactors = FALSE)
colnames(cluster_info) <- c("id", "cluster")

clinical.data <- read.delim(file.path(data_dir, "clinical.data.final.txt"),
                            stringsAsFactors = FALSE)
clinical.data$cluster <- NULL


##add cluster information to metadata
combined_metadata<-merge(huntington_metadata, sample.match, by.x="GEO_Accession (exp)", by.y="Sample_id")
metadata.for.cluster<- merge(combined_metadata, cluster_info, by.x = "Run", by.y= "id")
rownames(metadata.for.cluster)<-metadata.for.cluster$Run
#write.table(merged_clinical_data, "final-metadata.txt", row.names = T,quote = F,sep='\t')


##combine clinical data and metadata
merged_clinical_data<-merge(metadata.for.cluster, clinical.data, by.x="Run", by.y="Run")
rownames(merged_clinical_data)<-merged_clinical_data$Run




##export files
write.table(merged_clinical_data,
            file = file.path(out_dir, "cleaned-metadata.txt"),
            row.names = TRUE, quote = FALSE, sep = "\t")

write.csv(raw_count_data.filt,
          file = file.path(out_dir, "qc-count-data.csv"),
          row.names = TRUE, quote = FALSE)



