#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#argument 1 is path to scenic output
#argument 2 is path to IVF seurat object
#argument 3 is directory where files should be written out
#argument 4 is suffix to be added to files

# args=c("gc_scoring/auc_matrix_joined_firstbatch.csv",
#        "IVF.RData",
#        "gc_scoring/",
#        "first_batch")

library(Seurat)
library(ggplot2)
library(caret)
library(glmmTMB)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

auc_Train <- loadRData("auc_Train_new.Rdata") #created by script auc_classifier.R

#center and scaling
preProcValues <- preProcess(auc_Train, method = c("center", "scale"))

auc_df <- read.csv(args[1], 
                   row.names = 1)
auc_df <- as.data.frame(t(auc_df))
auc_df <- auc_df[,colnames(auc_df)%in%colnames(auc_Train)]
auc_df <- auc_df[grep("gc", rownames(auc_df)),]
auc_df <- auc_df[grep("so", rownames(auc_df)),]
auc_dfTransformed <- predict(preProcValues, auc_df)


############################
#gene matrix preparation

IVF <- loadRData(args[2])

seurat_m = IVF
seurat_m[["CellName"]] <- colnames(seurat_m)
exprs_matrix <- seurat_m@assays$SCT@data
exprs_matrix <- t(as.matrix(exprs_matrix))
exprs_matrix <- exprs_matrix[grep("so", rownames(exprs_matrix), ignore.case=T),]
exprs_matrix <- exprs_matrix[grep("gc", rownames(exprs_matrix), ignore.case=T),]
seurat_m <- subset(seurat_m, subset = CellName %in% rownames(exprs_matrix))
seurat_m <- SCTransform(seurat_m, return.only.var.genes = F)
exprs_matrix <- seurat_m@assays$SCT@data
exprs_matrix <- t(as.matrix(exprs_matrix))
exprs_matrix <- as.data.frame((exprs_matrix))

exprs_Train <- loadRData("exprs_Train_new.Rdata") #created by genes_classifier.R

# center and scaling
exprs_ivf <- exprs_matrix[grep("GC", rownames(exprs_matrix), invert=T),]
exprs_ivf <- exprs_ivf[,colnames(exprs_ivf)%in%colnames(exprs_Train)]
preProcValues <- preProcess(exprs_ivf, method = c("center", "scale"))
exprs_ivfTransformed <- predict(preProcValues, exprs_ivf)

setwd(args[3])

save(exprs_ivfTransformed, file=paste0("exprs_ivfTransformed", args[4], ".Rdata"))
save(auc_dfTransformed, file=paste0("auc_ivfTransformed", args[4], ".Rdata"))
##############################################################################
setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/IVF/classifier/final/")
load("auc_svm_2.Rdata") #model saved by auc_classifier.R
load("exprs_svm_new.Rdata") #model saved by genes_classifier.R

#predicting classes
svm_auc_ivf <- predict(auc_svm_2, newdata = auc_dfTransformed)
svm_exprs_ivf <- predict(exprs_svm, newdata = exprs_ivfTransformed)

#merging and saving
svm <- data.frame(cells=rownames(exprs_ivfTransformed),
                  auc=svm_auc_ivf,
                  exprs=svm_exprs_ivf)
setwd(args[3])
save(svm, file=paste0("gc_class", args[4], ".Rdata"))
