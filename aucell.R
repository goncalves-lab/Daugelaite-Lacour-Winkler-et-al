#this scripts calculates pathway activity scores

library(biomaRt)
library(AUCell)
library(stringr)
library(GSEABase)
library(easyGgplot2)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ggrepel)
#functions used to load gmt files
GeneSetListFromGmt <- function(path2gmt) {
  gene_sets <- lapply(readLines(path2gmt), str_split, "\t")
  set_list=lapply(gene_sets, handle_gene_sets) 
  return(set_list)
}

handle_gene_sets <- function(gene_sets) {
  gene_set_vec <- gene_sets[[1]]
  if (length(gene_set_vec[2]) == 0) 
    description <- "None" 
  else
    description <- gene_set_vec[2]
  genes <- gene_set_vec[3:length(gene_set_vec)][str_length(gene_set_vec[3:length(gene_set_vec)]) > 0]
  gene_set_list <- list(id = gene_set_vec[1],
                        description = description,
                        genes = genes)
  return( gene_set_list )
}

gene_set <- GeneSetListFromGmt("msigdb.v7.5.1.symbols.gmt")

pathways <- c()
for (i in seq(from=1, to=length(gene_set))){
  pathway_names <- gene_set[[i]][[1]]
  pathways <- c(pathway_names, pathways)
}

pathway_cc <- pathways[grep("CELL_CYCLE",pathways)]
pathway_sh <- pathways[grep("ANDROGEN",pathways)]
pathway_e <- pathways[grep("ESTROGEN",pathways)]
pathway_p <- pathways[grep("PROGESTERONE",pathways)]
pathway_ox <- pathways[grep("OXIDATIVE",pathways)]
pathway_an <- pathways[grep("FATTY",pathways)]
pathway_ch <- pathways[grep("CHROMOSOME",pathways)]
pathway_cp <- pathways[grep("CHECKPOINT",pathways)]
pathway_e2f <- pathways[grep("E2F",pathways)]
pathway_l <- pathways[grep("MAPK",pathways)]
pathway_g2m <- pathways[grep("G2_M",pathways)]
pathway_rep <- pathways[grep("REPAIR",pathways)]

pathway_selected <- data.frame(pathways=c(pathway_cc, pathway_sh, pathway_e, pathway_p, pathway_ox, pathway_an, 
                                 pathway_ch, pathway_cp, pathway_e2f, pathway_l, pathway_g2m, pathway_rep))

pathway_selected <- pathway_selected[c(66, 22, 107, 506, 139, 229, 575, 574, 608, 543, 122,39,16,182),]


pathways_list <- list()
for (i in seq(from=1, to=length(gene_set))){
  if (gene_set[[i]][[1]]%in%pathway_selected){
    pathway <- gene_set[i]
    pathways_list <- c(pathway, pathways_list) 
  }
  
}
#normalized data were produced usinf scnorm.R script
load("NormalizedData.Rdata")

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list <- getBM(attributes=c("external_gene_name", "ensembl_gene_id"), filters="ensembl_gene_id", values = rownames(NormalizedData), mart=ensembl)
NormalizedData <- merge(ensemnl_list, NormalizedData, by.x="ensembl_gene_id", by.y=0)
rownames(NormalizedData) <- make.unique(NormalizedData$external_gene_name)
NormalizedData <- NormalizedData[,3:193]
reads<-NormalizedData

rownames(reads) <- toupper(rownames(reads))
cells_rankings <- AUCell_buildRankings(as.matrix(reads), nCores=1, plotStats=TRUE)
geneSets_final <- list()
for (i in seq(from=1, to=length(pathways_list))){
  pathway=pathways_list[[i]][[3]]
  geneSets <- GeneSet(pathway, setName=pathways_list[[i]][[1]])
  geneSets_final <- c(geneSets, geneSets_final)
}


null_rank_final <- c()

for (i in seq(from=1, to=ncol(NormalizedData))){
  rank <- as.data.frame(NormalizedData[,i])
  #print(i)
  colnames(rank) <- "rank"
  rank <- rank[order(rank$rank, decreasing = T),]
  null_rank <- min(which(rank==0))
  #print(null_rank)
  null_rank_final <- c(null_rank_final, null_rank)
}


auc_df_final <- data.frame()
for(i in geneSets_final){
  print(i)
  cells_AUC <- try(AUCell_calcAUC(i, cells_rankings, aucMaxRank = ceiling(0.2 * nrow(cells_rankings))), silent=T)
  if (class(cells_AUC) != "try-error"){
    auc_df <- as.data.frame(getAUC(cells_AUC))
    #auc_df <- t(auc_df)
    auc_df_final <- rbind(auc_df, auc_df_final)}
}

auc_df_final <- t(auc_df_final)

metadata <- data.frame("condition"=colnames(reads))
metadata$condition <- gsub("_rep", "", metadata$condition)
metadata$condition <- substr(metadata$condition, 1, nchar(metadata$condition)-2)
metadata$condition <- gsub("S0", "SO", metadata$condition)
metadata$sample <- colnames(reads)
metadata$cell_type <- substr(metadata$condition, nchar(metadata$condition)-1, nchar(metadata$condition))
metadata$label <- substr(metadata$condition, 1, nchar(metadata$condition)-2)
metadata <- metadata[metadata$label%in%c("SO3M", "SO12M"),]

auc_df_final <- auc_df_final[rownames(auc_df_final)%in%metadata$sample,]
auc_df_final <- merge(metadata, auc_df_final, by.x="sample", by.y=0)

auc_df_final_oc <- auc_df_final[auc_df_final$cell_type=="GC",]

cells_all <-c()
cells_pathway <- list()
for(i in seq(from=5, to=ncol(auc_df_final_oc))){
  result_f <- c()
  auc_df <- auc_df_final_oc[,i]
  for(j in seq(length(auc_df))){
    result=sum(abs(auc_df) >= abs(auc_df[j]))/(length(auc_df))
    result_f <- c(result_f, result)
    
  }
  cells <- auc_df_final_oc[which(result_f<0.05),1]
  cells_all <- c(cells, cells_all)
  cells_pathway[[i-4]] <- cells 
}

names(cells_pathway) <- colnames(auc_df_final_oc)[5:ncol(auc_df_final_oc)]

auc_mean_final <- data.frame()
for (i in seq(from=5, ncol(auc_df_final_oc))){
  auc <- auc_df_final_oc[,c(4,i)]
  pathway_name <- colnames(auc)[2]
  colnames(auc)[2] <- "pathway"
  auc_mean <- auc %>% group_by(label) %>% summarize(Mean=mean(pathway))
  auc_mean$pathway <- pathway_name
  auc_mean_final <- rbind(auc_mean, auc_mean_final)
}


