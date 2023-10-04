library(Seurat)
library(ggplot2)
library(caret)
library(glmmTMB)

############################
#gene matrix preparation

genes <- read.table("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/cell_com/SO_quality_prediction/quality_marker_genes.tsv",
                    sep=",", header = T, row.names = 1)

load("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/data/20230523010_44/IVF.RData")
load("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/merged_19700_21475_22151/ovulation.Rdata")

seurat_m <- merge(IVF, ovulation, merge.data=F)
seurat_m[["CellName"]] <- colnames(seurat_m)
exprs_matrix <- seurat_m@assays$SCT@data
exprs_matrix <- t(as.matrix(exprs_matrix))
exprs_matrix <- exprs_matrix[grep("so", rownames(exprs_matrix), ignore.case=T),]
exprs_matrix <- exprs_matrix[grep("gc", rownames(exprs_matrix), ignore.case=T),]
exprs_matrix <- exprs_matrix[grep("12M", rownames(exprs_matrix), invert=T),]
seurat_m <- subset(seurat_m, subset = CellName %in% rownames(exprs_matrix))
seurat_m <- SCTransform(seurat_m)
exprs_matrix <- seurat_m@assays$SCT@data
exprs_matrix <- t(as.matrix(exprs_matrix))
exprs_matrix <- exprs_matrix[,colnames(exprs_matrix)%in%rownames(genes)]
exprs_matrix <- as.data.frame((exprs_matrix))



#identifying near-zero variance predictors
nzv <- nearZeroVar(exprs_matrix, saveMetrics= TRUE)

#identifying correlated predictors
descrCor <- cor(exprs_matrix)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .85)
exprs_matrix <- exprs_matrix[,-highlyCorDescr]

#identifying linear dependancies
comboInfo <- findLinearCombos(exprs_matrix)

#splitting dataset

load("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/merged_19700_21475_22151/clusters_Y_consensus.Rdata")
clusters_consensus <- clusters_consensus[!is.na(clusters_consensus)]
clusters_consensus <- clusters_consensus[grep("SO", names(clusters_consensus))]
clusters_consensus=droplevels(clusters_consensus)
clusters_df <- data.frame(cells=names(clusters_consensus), label=clusters_consensus)

exprs_so <- exprs_matrix[rownames(exprs_matrix)%in%clusters_df$cells,]
exprs_so <- exprs_so[match(clusters_df$cells, rownames(exprs_so)),]
set.seed(3456)
trainIndex <- createDataPartition(clusters_df$label, p = .7, 
                                  list = FALSE, 
                                  times = 1)
exprs_Train <- exprs_so[ trainIndex,]
exprs_Test  <- exprs_so[-trainIndex,]

#center and scaling
preProcValues <- preProcess(exprs_Train, method = c("center", "scale"))

exprs_trainTransformed <- predict(preProcValues, exprs_Train)
exprs_testTransformed <- predict(preProcValues, exprs_Test)

exprs_ivf <- exprs_matrix[grep("3M", rownames(exprs_matrix), invert=T),]
exprs_ivfTransformed <- predict(preProcValues, exprs_ivf)

setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/IVF/classifier")

save(exprs_ivfTransformed, file="exprs_ivfTransformed.Rdata")

exprs_trainTransformed <- merge(exprs_trainTransformed, clusters_df, by.x=0, by.y="cells")
rownames(exprs_trainTransformed) <- exprs_trainTransformed$Row.names
exprs_trainTransformed <- exprs_trainTransformed[,2:ncol(exprs_trainTransformed)]
save(exprs_trainTransformed, file="exprs_trainTransformed.Rdata")

exprs_testTransformed <- merge(exprs_testTransformed, clusters_df, by.x=0, by.y="cells")
rownames(exprs_testTransformed) <- exprs_testTransformed$Row.names
exprs_testTransformed <- exprs_testTransformed[,2:ncol(exprs_testTransformed)]
save(exprs_testTransformed, file="exprs_testTransformed.Rdata")
