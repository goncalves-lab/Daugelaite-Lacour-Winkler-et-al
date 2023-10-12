#Creating the Seurat object for the IVF experiment in mice

#-------Libraries---------

library(Seurat)


#-------Functions----------


#-------Data-----------

#Import data from ArrayExpress
metadata = read.csv(file = "metadata_filter.csv") #sample metadata table
reads_df = read.csv("count_table_all_samples.csv", row.names = 1)
reads_df = reads_df[, metadata$Name]

#-------Main---------

#Create Seurat object
IVF <- CreateSeuratObject(counts = reads_df, project = "IVF")
#Running normalisation and dimensionality reduction 
IVF <- SCTransform(IVF, verbose = FALSE,return.only.var.genes = FALSE)
IVF <- RunPCA(IVF, features = VariableFeatures(object = IVF))
IVF <- RunUMAP(IVF, dims = 1:15)
IVF <- FindNeighbors(IVF)
IVF <- FindClusters(IVF)

#Adding metadata to Seurat object
IVF$mouse = metadata$Individual
IVF$cell = metadata$Cell_type

save(IVF, file = "IVF.RData")


