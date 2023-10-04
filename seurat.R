#this script used Seurat to do QC filtering and normalize the data
library(Seurat)
library(biomaRt)
library(RColorBrewer)
library(ggrepel)
library(dplyr)

color_vector <- c("#31135e","#6a29c4","#be95ff","#e8daff","#510224","#9f1853",
                  "#ff7eb6","#ffd6e8")
file_readeR <- function(file_dir){
  file_names <- list.files(file_dir)
  reads_log <- read.table(file_names[1], sep = "\t",fill = T)
  reads_log <- reads_log[5:nrow(reads_log),1:2]
  colnames(reads_log) <- c("gene_names", gsub("_ReadsPerGene.out.tab","",file_names[1]))
  for (i in seq(from=2, to=length(file_names))){
    reads_log_temp <- read.table(file_names[i], sep = "\t",fill = T)
    reads_log_temp <- as.data.frame(reads_log_temp[5:nrow(reads_log_temp),2])
    colnames(reads_log_temp) <- c(gsub("_ReadsPerGene.out.tab","",file_names[i]))
    reads_log <- cbind (reads_log, reads_log_temp)
  }
  return(reads_log)
}


reads_df <- file_readeR("gene_counts/")
reads_df$gene_names <- gsub("\\..*","", reads_df$gene_names)
rownames(reads_df) <- reads_df$gene_names

sample_list <- read.csv("metadata.tsv", header=T, sep="\t")
sample_list <- sample_list[,c(1,11)]
sample_list$FASTQ_FILE <- gsub("_R1.fastq.gz","",sample_list$FASTQ_FILE)
colnames(sample_list) <- c("fastq", "sample")

reads_df=as.data.frame(t(reads_df))
reads_df$sample_names <- rownames(reads_df)
reads_df <- merge(reads_df, sample_list, by.x = "sample_names", by.y = "fastq" )
rownames(reads_df) <- reads_df$sample

reads_df <- as.data.frame(t(reads_df))
reads_df <- reads_df[-c(1,55489),]  

qc <- read.csv("metadata/qc_table.csv", sep=",", header=T)
qc_genes <- read.csv("metadata/gene_qc_table.csv", sep=",", header = T)

qc <- qc[qc$Number_of_reads>5000,]
qc_metrics <- merge(qc, qc_genes, by.x="samples",by.y = "sample")
qc_metrics <- qc_metrics[qc_metrics$per_mt_rna<5,]
qc_metrics <- qc_metrics[qc_metrics$number_of_genes>1000,]

sample_vec <- qc_metrics$samples
reads_df <- reads_df[,sample_vec]
reads_df<- reads_df[colnames(reads_df)%in%qc_metrics$samples,]

ovulation <- CreateSeuratObject(counts = reads_df, project = "ovulation")

ovulation <- SCTransform(ovulation, verbose = FALSE,return.only.var.genes = FALSE)
ovulation <- RunPCA(ovulation, features = VariableFeatures(object = ovulation))
ovulation <- RunUMAP(ovulation, dims = 1:30)
ovulation <- FindNeighbors(ovulation, dims = 1:30)
ovulation <- FindClusters(ovulation, resolution = 0.5)

