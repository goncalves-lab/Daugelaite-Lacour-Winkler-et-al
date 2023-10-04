#normalization of data using SCnorm package
library(SCnorm)
library(SingleCellExperiment)
library(dplyr)
library(biomaRt)
library(Seurat)

#reads_df <- read.csv("reads_all_samples.csv", row.names = 1)
load("ovulation.Rdata") #produced in script seurat.R
reads_df <- as.data.frame(ovulation[["RNA"]]@counts)

metadata <- data.frame("condition"=colnames(reads_df))
metadata$condition <- gsub("_rep", "", metadata$condition)
metadata$condition <- substr(metadata$condition, 1, nchar(metadata$condition)-2)
metadata$condition <- gsub("S0", "SO", metadata$condition)
metadata <- metadata %>% mutate(label = case_when(
  condition == "SO3MOC" ~ 1,
  condition == "NO3MOC" ~ 2,
  condition == "SO12MOC" ~ 3,
  condition == "NO12MOC" ~ 4,
  condition == "SO3MGC" ~ 5,
  condition == "NO3MGC" ~ 6,
  condition == "SO12MGC" ~ 7,
  condition == "NO12MGC" ~ 8))
metadata$sample <- colnames(reads_df)

metadata$counts <- colSums(reads_df)
metadata <- metadata[metadata$counts>200,] 
reads_df <- reads_df[,colnames(reads_df)%in%metadata$sample]  


gene_info <- data.frame("ensembl_id"=rownames(reads_df))

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
filters <- listFilters(ensembl)
atributes <- listAttributes(ensembl) 
ensemnl_list <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "transcript_length", "transcript_appris"), filters="ensembl_gene_id", values = gene_info$ensembl_id, mart=ensembl)

ensemnl_list <- ensemnl_list[grep("principal", ensemnl_list$transcript_appris),]

ensemnl_list <- ensemnl_list %>%
  group_by(ensembl_gene_id) %>%
  summarize(Mean = mean(transcript_length, na.rm=TRUE))

genes_missing <- as.data.frame(gene_info[!gene_info$ensembl_id%in%ensemnl_list$ensembl_gene_id,])
ensemnl_list_missing <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "transcript_length", "transcript_appris"), filters="ensembl_gene_id", values = genes_missing$`gene_info[!gene_info$ensembl_id %in% ensemnl_list$ensembl_gene_id, ]`, mart=ensembl)

ensemnl_list_missing <- ensemnl_list_missing %>%
  group_by(ensembl_gene_id) %>%
  summarize(Mean = mean(transcript_length, na.rm=TRUE))

ensemnl_list <- rbind(ensemnl_list, ensemnl_list_missing)

full_list <- merge(ensemnl_list, reads_df,by.x="ensembl_gene_id", by.y=0)
gene_length <- full_list$Mean

rownames(full_list) <- full_list$ensembl_gene_id


reads_df <- full_list[,3:193]
ovulation_sce <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = reads_df))

Conditions <- metadata$label


countDeptEst <- plotCountDepth(Data = ovulation_sce, Conditions = Conditions,
                               FilterCellProportion = .1, NCores=3)



DataNorm <- SCnorm(Data = ovulation_sce,
                   Conditions = Conditions,
                   PrintProgressPlots = TRUE,
                   FilterCellNum = 10,
                   NCores=1, reportSF = TRUE,
                   withinSample = gene_length)


NormalizedData <- SingleCellExperiment::normcounts(DataNorm)
save(NormalizedData, file="NormalizedData_expanded.Rdata")


