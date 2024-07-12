#Pseudotime analysis for oocytes and granulosa cells
#This script is not used anymore for the final version of the manuscript


#--------Libraries---------

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(uwot)
library(mclust, quietly = TRUE)
library(RColorBrewer)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))

#-------Functions----------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

ensembleR=function(val,attr,filt){
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
  filters=listFilters(ensembl)
  attributes=listAttributes(ensembl)
  ensembl_df=getBM(attributes=attr,filters=filt,values = val,mart=ensembl)
  return(ensembl_df)
}

subseteR<- function(seurat_obj, gene_subset){
  cell_counts <-  as.matrix(seurat_obj[["SCT"]]@counts)
  cell_counts <- cell_counts[rownames(cell_counts)%in%gene_subset,]
  
  colnames(cell_counts) <- gsub("_rep", "", colnames(cell_counts))
  colnames(cell_counts) <- gsub("S012M", "SO12M", colnames(cell_counts))
  sce <- SingleCellExperiment(list(logcounts = cell_counts))
  return(sce)
}

reductoR <- function(sce){ #runs different dimensionality reductions and clusterings
  pca <- prcomp(t(log1p(assays(sce)$logcounts)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  print(plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1))
  rd2 <- umap(t(log1p(assays(sce)$logcounts)))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  print(plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1))
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  print(plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1))                               
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(sce)$kmeans <- cl2
  print(plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1))
  return(sce)
}

label_makeR <- function(sce){
  label_df <- as.data.frame(rownames(reducedDims(sce)$PCA))
  label_df$label <- substr(label_df$`rownames(reducedDims(sce)$PCA)`,1,nchar(label_df$`rownames(reducedDims(sce)$PCA)`)-2)
  
  label_df %<>% mutate(color_group = case_when(
    label == "NO12MOC" ~ brewer.pal(11, "Spectral")[10],
    label == "NO3MOC" ~ brewer.pal(11, "Spectral")[8],
    label == "SO12MOC" ~ brewer.pal(11, "Spectral")[2],
    label == "SO3MOC" ~ brewer.pal(11, "Spectral")[5]
  ))
  return(label_df)
}


ploteR <- function(sce){ #plot the results from slingshot
  curve_sce <- as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
  pca_pt <- as.data.frame(reducedDims(sce)$PCA)
  pca_pt$pseudotime_1 <-  sce$slingPseudotime_1
  print(ggplot(pca_pt, aes(x=PC1, y=PC2, color= pca_pt$pseudotime_1)) + geom_point(size=5)+
          theme_bw(base_size = 14) + 
          theme(legend.position = "right")+
          theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x="PC1",y="PC2")+
          theme(axis.text=element_text(size=34), axis.title=element_text(size=36))+
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+
          geom_smooth(aes(x=PC1, y=PC2), data=curve_sce, level=0, size=2)+
          scale_colour_gradientn(colors = brewer.pal(11, "Spectral")[11:1]))
  pca_pt$label <- row.names(pca_pt)  
  pca_pt$group_label <- substr(pca_pt$label, 1, nchar(pca_pt$label)-2)
  print(ggplot(pca_pt, aes(x=PC1, y=PC2, color= pca_pt$group_label)) + geom_point(size=5)+
          theme_bw(base_size = 14) + 
          theme(legend.position = "right")+
          theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x="PC1",y="PC2")+
          theme(axis.text=element_text(size=34), axis.title=element_text(size=36))+
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+
          geom_smooth(aes(x=PC1, y=PC2), data=curve_sce, level=0, color="black", size=2)+
          scale_colour_manual(values=brewer.pal(11, "Spectral")[c(2,4,8,11)]))
  print(ggplot(pca_pt, aes(x= group_label, y=pseudotime_1)) + geom_point(size=5, color="black")+
          theme_bw(base_size = 14) + 
          theme(legend.position = "right")+
          theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x="", y="pseudotime")+
          theme(axis.text=element_text(size=26, angle = 0), axis.title=element_text(size=28))+
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+ coord_flip()+
          scale_colour_manual(values=brewer.pal(11, "Spectral")[c(2,4,8,11)]))
  return(pca_pt)
}

ks_testeR <- function(pt, comp_list){ #test difference between groups
  res <- data.frame()
  for (i in seq(from=1, to=length(comp_list))){
    condition_a <- pt[pt$group_label == comp_list[[i]][1],]
    condition_b <- pt[pt$group_label == comp_list[[i]][2],]
    dst <- ks.test(condition_a$pseudotime_1, condition_b$pseudotime_1)
    res_temp <- data.frame("con_1" = comp_list[[i]][1], "con_2" = comp_list[[i]][2], "p_value" = dst[[2]])
    res <- rbind(res, res_temp)
  }
  res$p_adj <- p.adjust(res$p_value, method = "BH", n = length(res$p_value))  
  return(res)
}

#--------Data---------

ovulation<- loadRData("ovulation.Rdata") #Seurat object created by script create_seurat_age_ov.R
Idents(ovulation) <- ovulation$cell
ovulation_oc <- subset(ovulation, idents="GC") #Here select the cell type you want to study (OC or GC)

#-------Main----------

var_genes <- ovulation_oc@assays$SCT@var.features

sce_hvg <- subseteR(ovulation_oc, var_genes)

sce_hvg <- reductoR(sce_hvg)

sce_hvg <- slingshot(sce_hvg, clusterLabels = 'GMM', reducedDim = 'PCA')

pca_pt <- as.data.frame(reducedDims(sce_hvg)$PCA)

pt_hvg <- ploteR(sce_hvg)


comp_list <- list(c("NO3MOC","NO12MOC"), c("NO3MOC", "SO3MOC"), c("NO12MOC", "SO12MOC"),c("SO3MOC", "SO12MOC"))
res_hvg <- ks_testeR(pt_hvg, comp_list)


