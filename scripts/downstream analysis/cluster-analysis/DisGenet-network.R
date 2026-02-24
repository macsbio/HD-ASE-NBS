library(RCy3)
library(httr)
#install.packages("devtools")
library(devtools)
#install_gitlab("medbio/disgenet2r")
library(disgenet2r)
library(biomaRt)

# set working dir to the folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# paths relative to the script directory
data_dir <- file.path(getwd(), "Data")
out_dir  <- file.path(getwd(), "output")

# create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# (optional) stop if Data folder is missing
if (!dir.exists(data_dir)) stop("Missing Data folder: ", data_dir)

key.genes <- read.delim(file.path(out_dir, "sig.genes.ase.txt"), sep ="\t", header = F)


library(disgenet2r)
api_key <- "351515e2-f987-4f31-ab92-1c9fa3f5a534"
Sys.setenv(DISGENET_API_KEY= api_key)

hd <- disease2gene(
  disease   = "UMLS_C0020179",   # Huntington’s disease (UMLS CUI)
  database  = "CURATED",    # with academic license you must stick to CURATED
  score     = c(0, 1)
)

hd_genes <- data.frame(gene=unique(hd@qresult[["gene_symbol"]]),
                      score = hd@qresult[["score"]],
                      evidence_index = hd@qresult[["evidence_index"]]) 



#convert to gene symbols
grch37.gene <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")


##gene categories 
genes<- getBM(mart=grch37.gene, values= hd_genes$gene,
              filters= "hgnc_symbol",
              attributes= c("ensembl_gene_id", "hgnc_symbol",
                            "gene_biotype"))

combined_gene_dataframe<-data.frame(gene = c(genes$ensembl_gene_id, key.genes$V1),
                                    groups = c(replicate(nrow(genes), "HD genes"),replicate(nrow(key.genes), "key genes")))


write.table(combined_gene_dataframe$gene, "Disgenet-network.txt", quote = F, row.names = F, col.names = F)
write.table(combined_gene_dataframe, "Disgenet-table-network.txt", quote = F, row.names = F, col.names = T, sep = " ")


# 1. getSTRINGdb for human
# Apply 0.9 confidence score
string_db_all <- STRINGdb$new(
  version         = "12.0",  # or "11"
  species         = 9606,
  score_threshold = 0
)

# Get the whole network as an igraph object
g_all <- string_db_all$get_graph()

# Convert to an edge list data.frame
ppi_all <- igraph::as_data_frame(g_all, what = "edges")
head(ppi_all)
dim(ppi_all)

# 3. map gene ids to protein ids
### extract protein ids from the human network
protein_ids <- sapply(strsplit(ppi_all$from, '\\.'),
                      function(x) x[2])

### get protein to gene id mappings
mart_results <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                      filters = "ensembl_peptide_id", 
                      values = protein_ids,
                      mart = grch37.gene)
### extract protein ids from the human network
ppi_all$from <- sapply(strsplit(ppi_all$from, '\\.'),
                      function(x) x[2])
### extract protein ids from the human network
ppi_all$to <- sapply(strsplit(ppi_all$to, '\\.'),
                      function(x) x[2])

## 1. Make sure mart_results has one row per peptide id
mart_clean <- mart_results %>%
  distinct(ensembl_peptide_id, .keep_all = TRUE)

protein_id_match<- match(ppi_all$from, mart_clean$ensembl_peptide_id)
ppi_all$from_gene_id <- mart_clean$ensembl_gene_id[protein_id_match]
protein_id_match_to<- match(ppi_all$to, mart_clean$ensembl_peptide_id)
ppi_all$to_gene_id <- mart_clean$ensembl_gene_id[protein_id_match_to]


##remove NAs
ppi_all<- ppi_all[complete.cases(ppi_all),]

ppi_all_filtered <- ppi_all %>%
  filter(from_gene_id %in% combined_gene_dataframe$gene)




## HD gene set
hd_ids <- merged_hd_genes$ensembl_gene_id

## 1) Edges where *at least one* endpoint is an HD gene
ppi_hd_any <- ppi_with_genes %>%
  filter(gene_from %in% hd_ids | gene_to %in% hd_ids)

##subset those interactions ith score of above 400
ppi_hd_any<- subset(ppi_hd_any, combined_score >= 900)

cytoscapePing()

gene_nodes <- data.frame(id=pathway_id$id,
                                 group= pathway_id$group,
                                 stringsAsFactors = FALSE)

loadTableData(combined_gene_dataframe, data.key.column = "gene", table.key.column = "Ensembl")



# Generate the gene-disease network
geneDisResult <- disgenetRestCall("gene-disease-net",geneDisParams)

layoutNetwork("force-directed", network = as.numeric(geneDisResult$result$networkSUID))

title <- "mentaldisorder-GDA-network"
renameNetwork(title, network =   as.numeric(geneDisResult$result$networkSUID))

huntington_genes_disg<-getTableColumns()
huntington_genes_disg<-huntington_genes_disg[-c(1),]

##combine genes ( key and huntington genes)

combined.key.genes<- data.frame(genes = c(genes$hgnc_symbol, disgenet_data$Gene),
                                group = c(replicate(nrow(genes), "Key genes"),replicate(nrow(disgenet_data), "HD-associated genes")))

loadTableData(combined.key.genes)
loadTableData(
  combined.key.genes,
  data.key.column = "query_term",
  table = "node",
  table.key.column = "query term",
  namespace = "default",
  network = "current")



disgenet_data_curated<- read.delim("DISEASES_Summary_GDA_CURATED_C0020179.tsv")
disgenet_data_all<- read.delim("DISEASES_Summary_GDA_ALL_C0020179.tsv")
