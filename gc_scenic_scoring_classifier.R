#!/usr/bin/env Rscript
#Computing AUC score based on reference SCENIC regulons

args = commandArgs(trailingOnly=TRUE)

#argument 1 is path to count table normalised with SCnorm
#make sure that gc samples are named in this manner: e.g. so75gc1, otherwise line 93 
#has to be modified
#argument 2 is directory where files should be written out
#argument 3 is suffix to be added to files

# args=c("NormalizedData_expanded.Rdata", "gc_scoring/", "")

library(SCENIC)
library(gsubfn) 
library(biomaRt)
library(Seurat)
library(ggplot2)
library(e1071)
library(glmmTMB)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

########################################################

input_file = args[[1]]
output_dir = args[[2]]
suffix = args[[3]]

# loading data of samples to be classified, this data should be a count table normalised by Scnorm.R
ovulation = loadRData(input_file)

exprs_matrix <- ovulation

#selecting granulosa cells
exprs_matrix <- exprs_matrix[,grep("gc", colnames(exprs_matrix), ignore.case=T)]

#Changing gene names
cell_info <- data.frame(label = colnames(exprs_matrix))
cell_info$nGene <- colSums(exprs_matrix >0)
rownames(cell_info) <- make.names(cell_info$label, unique = T)
ensembl=suppressWarnings(useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl"))
ensembl_list=getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values = row.names(exprs_matrix), mart=ensembl)

for(i in seq(from=1, to=nrow(ensembl_list))){
  if (ensembl_list[i,2]==""){
    ensembl_list[i,2] = "no_id"
  }
}
full_list=merge(exprs_matrix, ensembl_list,by.x=0,by.y="ensembl_gene_id")
rownames(full_list) <- make.names(full_list$external_gene_name, unique=T)
exprs_matrix <- as.matrix(full_list[,2:(ncol(exprs_matrix)+1)])

cell_info <- cell_info[grep("gc", cell_info$label, ignore.case = T),]
exprs_matrix <- exprs_matrix[,colnames(exprs_matrix)%in%cell_info$label]


#Selecting the same genes that were used by Scenic on the reference dataset (normalised with SCnorm)
genesKept = readRDS("int/1.1_genesKept.Rds")
exprMat_filtered <- exprs_matrix[rownames(exprs_matrix)%in%genesKept, ]
exprMat_filtered <- log2(exprMat_filtered+1)

#Scoring the regulons that were used by Scenic on reference dataset (normalised with SCnorm)
regulons <- readRDS("int/3.1_regulons_forAUCell.Rds")
aucellRankings = AUCell::AUCell_buildRankings(exprMat_filtered, plotStats=TRUE)
regulonAUC = AUCell::AUCell_calcAUC(regulons, aucellRankings, 
                                    aucMaxRank=2000)

auc_matrix <- AUCell::getAUC(regulonAUC)

#filtering the results based on reference Scenic run (normalised with SCnorm)
#this file is obtained by running this script on the reference dataset and skipping these two lines
auc_ref = read.csv("auc_matrix_joined_ref.csv", row.names = 1) 
auc_matrix = auc_matrix[rownames(auc_matrix)%in%rownames(auc_ref),]

setwd(output_dir)

save(exprs_matrix, file=paste0("exprs_matrix", suffix, ".Rdata"))
save(cell_info, file=paste0("cell_info", suffix, ".Rdata"))

write.csv(auc_matrix, paste0("auc_matrix_joined", suffix, ".csv"), quote = F)



