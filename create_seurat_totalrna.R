#Create total-RNA seq Seurat object

#--------Libraries---------

library(Seurat)
library(biomaRt)

#--------Data--------------

#Import data from ArrayExpress
count_table = read.csv("count_table_all_samples.csv", row.names = 1)
metadata = read.csv("sample_table.csv") #sample metadata table

#--------Main--------------

#Create Seurat object and add metadata
tot_RNA = CreateSeuratObject(count_table)
tot_RNA$cell = c("OC", "GC")[1+grepl("GC", colnames(tot_RNA))]
tot_RNA$ov =  substr(colnames(tot_RNA), 1,2)
tot_RNA$mouse = metadata$Individual[match(colnames(tot_RNA), metadata$Name)]
#Run normalisation and dimensonality reductions
tot_RNA = SCTransform(tot_RNA, return.only.var.genes = F)
tot_RNA = RunPCA(tot_RNA, npcs = 10)
tot_RNA = FindNeighbors(tot_RNA, dims = 1:5)
tot_RNA = FindClusters(tot_RNA)
tot_RNA = RunUMAP(tot_RNA, dims = 1:5)

save(tot_RNA, file = "tot_RNA.RData")
