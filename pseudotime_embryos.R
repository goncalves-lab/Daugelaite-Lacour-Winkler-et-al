#Pseudotime analysis of mouse embryos



#-------Libraries-----------

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(uwot)
library(mclust, quietly = TRUE)
library(biomaRt)
suppressPackageStartupMessages(library(magrittr))
library(pheatmap)
library(factoextra)

#-------Functions-----------

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
  sce <- SingleCellExperiment(list(logcounts = cell_counts))
  return(sce)
}

reductoR <- function(sce){
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
  colData(sce)$group = as.factor(embryos$group[colnames(sce)])
  print(plot(rd1, col = brewer.pal(9,"Set1")[as.factor(embryos$group[colnames(sce)])], pch=16, asp = 1))
  return(sce)
}


ploteR <- function(sce){
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
  
  pca_pt$group_label = embryos$group[rownames(pca_pt)]
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
  
  print(ggplot(pca_pt, aes(x= group_label, y=pseudotime_1, color= group_label)) + geom_point(size=5)+
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

ks_testeR <- function(pt, comp_list){
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


#-------Data------------

IVF<- loadRData("IVF.RData") #Seurat object created by script XXXXXX
Idents(IVF) <- IVF$cell
embryos <- subset(IVF, idents="morula")
embryos = SCTransform(embryos)
embryos = RunPCA(embryos, npcs = 10)
stopped = c("so75e3", "so75e6", "so75e12", "so75e17", "so76e14", "so80e14", "so80e20", "so81e4")
embryos$development = ifelse(colnames(embryos)%in%stopped, "stopped", "morula")

load("gc_class.Rdata") # class prediction from script XXXX
names(gc_class) = gsub("gc", "e", names(gc_class))
embryos$group = gc_class[colnames(embryos)]
embryos$group[embryos$development=="stopped"] = "stopped"

#-------Main-----------

ensembl_mouse=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list_mouse=getBM(attributes=c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=ensembl_mouse)


#variable genes
var_genes <- embryos@assays$SCT@var.features

sce_hvg <- subseteR(embryos, var_genes)

sce_hvg <- reductoR(sce_hvg)

sce_hvg <- slingshot(sce_hvg, clusterLabels = 'GMM', reducedDim = 'PCA')

pt_hvg <- ploteR(sce_hvg)


to_rm = which(pt_hvg$pseudotime_1<20 & pt_hvg$group_label=="YS_good")
stripchart(pt_hvg$pseudotime_1[-to_rm] ~ factor(pt_hvg$group_label[-to_rm], levels = c("YS_good", "NA", "YS_bad", "stopped"), labels = c("S1", "unclassified", "S2", "stopped\ndeveloping")),
           pch = 20, method = "jitter", xlab = "pseudotime", las = 1)
segments(x0 = aggregate(pseudotime_1 ~ group_label, pt_hvg[-to_rm,], mean)$pseudotime_1, y0 = c(1.8, 3.8, 2.8, 0.8), y1 = c(2.2, 4.2, 3.2, 1.2), col = "red")

comp_list <- list(c("YS_bad","YS_good"))
res_hvg <- ks_testeR(pt_hvg, comp_list)



