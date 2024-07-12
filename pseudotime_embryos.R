#Pseudotime analysis of mouse embryos



#-------Libraries-----------

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(uwot)
library(biomaRt)
library(mclust, quietly = TRUE)
suppressPackageStartupMessages(library(magrittr))
library(RColorBrewer)

#-------Functions-----------

subseteR<- function(seurat_obj, gene_subset){
  cell_counts <-  as.matrix(seurat_obj[["SCT"]]@counts)
  cell_counts <- cell_counts[rownames(cell_counts)%in%gene_subset,]
  sce <- SingleCellExperiment(list(logcounts = cell_counts))
  return(sce)
}

reductoR <- function(sce, seuratobj, grouping){ #runs different dimensionality reductions and clusterings
  pca <- prcomp(t(log1p(assays(sce)$logcounts)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  print(plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1))
  reducedDims(sce) <- SimpleList(PCA = rd1)
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  print(plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1))                               
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(sce)$kmeans <- cl2
  print(plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1))
  colData(sce)$stage = as.factor(seuratobj@meta.data[colnames(sce), grouping])
  print(plot(rd1, col = brewer.pal(9,"Set1")[as.factor(seuratobj@meta.data[colnames(sce), grouping])], pch=16, asp = 1))
  return(sce)
}


ploteR <- function(sce, out_str, seuratobj, grouping){ #plot the results from slingshot
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
  
  pca_pt$group_label = seuratobj@meta.data[rownames(pca_pt), grouping]
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
          scale_colour_manual(values=brewer.pal(12, "Spectral")[c(2,4,6,8,10, 11,12)]))
  
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
          scale_colour_manual(values=brewer.pal(12, "Spectral")[c(2,4,6,8,10, 11,12)]))
  return(pca_pt)
}

dif_testeR <- function(pt, comp_list){ #test difference between groups
  res <- data.frame()
  for (i in seq(from=1, to=length(comp_list))){
    condition_a <- pt[pt$group_label == comp_list[[i]][1],]
    condition_b <- pt[pt$group_label == comp_list[[i]][2],]
    dst <- wilcox.test(condition_a$pseudotime_1, condition_b$pseudotime_1)
    res_temp <- data.frame("con_1" = comp_list[[i]][1], "con_2" = comp_list[[i]][2], "p_value" = dst$p.value)
    res <- rbind(res, res_temp)
  }
  res$p_adj <- p.adjust(res$p_value, method = "BH", n = length(res$p_value))  
  return(res)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}



#-------Data------------

ensembl_mouse=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list_mouse=getBM(attributes=c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=ensembl_mouse)


# Importing reference dataset from Xue et al 2013 (obtained from ArrayExpress)
reads_df = read.csv("reads_all_samples.csv", row.names = 1)
metadata = read.csv("metadata_filter.csv")
metadata$Sample_name = gsub(" ", ".", metadata$Sample_name)
metadata$Sample_name = gsub("-", ".", metadata$Sample_name)



# Importing experimental IVF data
IVF = loadRData("IVF.RData") #Seurat object created by script create_seurat_ivf_mouse.R
Idents(IVF) = IVF$cell
embryos = subset(IVF, idents!="granulosa")
embryos = SCTransform(embryos)
embryos = RunPCA(embryos, npcs = 10)

embryo_groups = read.csv("metadata.csv", header = T)[, c("Name", "Developmental Stage")] # metadata table from ArrayExpress
stages = c("so75" = "morula", "so76" = "morula", "so80" = "morula", "so81" = "morula", 
           "so88" = "blastocyst", "so90" = "blastocyst", 
           "so91" = "blastocyst", "so92" = "blastocyst", 
           "no34" = "blastocyst", "no35" = "blastocyst", "no36" = "blastocyst", 
           "no37" = "blastocyst", "no39" = "blastocyst", "no40" = "blastocyst", 
           "no38" = "morula", "no41" = "morula", "no42" = "morula")
embryo_groups$stage = stages[substr(embryo_groups$Name)]
colnames(embryo_groups) = c("embryo", "dev", "stage")
rownames(embryo_groups) = embryo_groups$embryo



load("gc_class.Rdata") # class prediction from classifier
names(gc_class) = gsub("gc", "e", names(gc_class))
embryos$group = gc_class[colnames(embryos)]
embryos$group[embryos$development=="stopped"] = "stopped"
embryos$group[embryos$group=="NA"] = "unclassified"
embryos$group[embryos$development=="unclear"] = NA

#-------Main-----------

## Generating the reference developmental dataset based on Xue et al data and our IVF embryos derived from natural ovulation
nat_counts = as.data.frame(IVF@assays$RNA@counts[,IVF$ovulation == "no" & IVF$cell != "granulosa" & embryo_groups[colnames(IVF), "dev"]%in%c("morula", "blastocyst")])
reads_df = merge(reads_df, nat_counts, by = 0)
rownames(reads_df) = reads_df$Row.names
reads_df = reads_df[,-1]
metadata = rbind(metadata, data.frame(Sample_name = rownames(IVF@meta.data[IVF$ovulation == "no" & IVF$cell != "granulosa" & embryo_groups[colnames(IVF), "dev"]%in%c("morula", "blastocyst"),]),
                                      Stage = IVF$cell[IVF$ovulation == "no" & IVF$cell != "granulosa" & embryo_groups[colnames(IVF), "dev"]%in%c("morula", "blastocyst")],
                                      Sample_number = NA, ENA_exp = NA, ENA_run = NA))

timecourse <- CreateSeuratObject(counts = reads_df, project = "reference_Huang2013")

timecourse <- SCTransform(timecourse, verbose = FALSE,return.only.var.genes = FALSE)
timecourse <- RunPCA(timecourse, features = VariableFeatures(object = timecourse), npcs=5)

timecourse$stage = metadata$Stage[match(names(timecourse$orig.ident), metadata$Sample_name)]
timecourse$mitopercent = qc$per_mt_rna[match(names(timecourse$orig.ident), qc$sample)]
timecourse$replicate = ifelse(grepl("no", colnames(timecourse)), "nat IVF", "Xue2013")
DimPlot(timecourse, group.by = "stage", pt.size = 5)
DimPlot(timecourse, group.by = "replicate", pt.size = 5)


var_genes <- timecourse@assays$SCT@var.features


## Finding genes that are linked to developmental stage


colData = data.frame(row.names = names(timecourse$orig.ident),
                     stage = timecourse$stage,
                     score = c("oocyte" = 1, "pronucleus" = 2, "2-cell blastomere" = 3, # stages are ordered
                               "4-cell blastomere" = 4, "8-cell blastomere" = 5, "morula" = 6, "blastocyst" = 7)[timecourse$stage])

#getting the variable genes expression
expr_table = timecourse@assays$SCT@counts[var_genes,]
expr_table = expr_table[rowMeans(expr_table)>1,]

#merging the expression table with the stage information
model_table = cbind(colData, t(expr_table[,]))
model_table = model_table[!model_table$stage%in%c("oocyte", "pronucleus"),] #excluding stages before ZGA

#computing Spearman correlation for each gene between its expression and the stage
cor_table = as.data.frame(t(apply(model_table[,3:ncol(model_table)], 2, 
                                  function(x){unlist(cor.test(x =x, y = model_table$score, method = "spearman")[c("estimate", "p.value")])})))
cor_table$gene_name = ensemnl_list_mouse$external_gene_name[match(rownames(cor_table), ensemnl_list_mouse$ensembl_gene_id)]
cor_table$padj = p.adjust(cor_table$p.value, method = "BH") #multiple testing correction

#filtering the genes based on significance and correlation, and on strong change between 8-cell - morula - blastocyst stages
dev_genes = na.exclude(rownames(cor_table)[cor_table$padj<0.01 & abs(cor_table$estimate.rho)>0.85])
fc_8tomorula = log2(colMeans(mlm_table[mlm_table$stage=="morula", 3:ncol(mlm_table)]) / 
                      colMeans(mlm_table[mlm_table$stage=="8-cell blastomere", 3:ncol(mlm_table)]))
fc_morutoblasto = log2(colMeans(mlm_table[mlm_table$stage=="blastocyst", 3:ncol(mlm_table)]) / 
                         colMeans(mlm_table[mlm_table$stage=="morula", 3:ncol(mlm_table)]))
dev_genes = dev_genes[sapply(dev_genes,
                             function(x){
                               return(abs(fc_8tomorula[x])>1 & abs(fc_morutoblasto[x])>1 & 
                                        fc_8tomorula[x] * cor_table[x,"estimate.rho"] > 0 & fc_morutoblasto[x] * cor_table[x,"estimate.rho"] > 0)
                             })]
dev_genes = na.exclude(dev_genes)

write.csv(cor_table[dev_genes, ], file = "developmental_pseudotime_genes.csv", quote = F)


## Validation of the selected genes to predict embryonic developmental stage

#running the pseudotime analysis on the reference dataset using the selected genes
sce_hvg_val <- subseteR(timecourse, dev_genes)
sce_hvg_val <- reductoR(sce_hvg_val, timecourse, "stage")
sce_hvg_val <- slingshot(sce_hvg_val, reducedDim = 'PCA')
pt_hvg_val <- ploteR(sce_hvg_val, "oc_hvg", timecourse, "stage")



#computing correlation as before for comparison
colData = data.frame(row.names = names(timecourse$orig.ident),
                     stage = timecourse$stage,
                     score = pt_hvg_val$pseudotime_1[match(names(timecourse$orig.ident), rownames(pt_hvg_val))])

expr_table = timecourse@assays$SCT@counts[dev_genes,]
expr_table = expr_table[rowMeans(expr_table)>1,]

mlm_table_val = cbind(colData, t(expr_table[,]))
mlm_table_val = mlm_table_val[!mlm_table_val$stage%in%c("oocyte", "pronucleus"),]

cor_table_val = as.data.frame(t(apply(mlm_table_val[,3:ncol(mlm_table_val)], 2, function(x){unlist(cor.test(x =x, y = mlm_table_val$score, method = "spearman")[c("estimate", "p.value")])})))
cor_table_val$gene_name = ensemnl_list_mouse$external_gene_name[match(rownames(cor_table_val), ensemnl_list_mouse$ensembl_gene_id)]
cor_table_val$padj = p.adjust(cor_table_val$p.value, method = "BH")



## Running pseudotime analysis on all IVF embryos

embryos$ovulation = substr(colnames(embryos), 1,2)
embryos$development = embryo_groups[colnames(embryos), "dev"]


#Adding information in the Seurat object about the groups we want to compare
#in the pseudotime, i.e. natural, SN, S, and arrested
embryos$group = gc_class[colnames(embryos)]
embryos$group[embryos$ovulation=="no"] = "nat"
embryos$group[embryos$development%in%c("early", "middle", "late", "oocyte")] = "stopped"
embryos$group[embryos$development=="unclear"] = "unclear"

DimPlot(embryos, group.by = "group", reduction = "pca", pt.size = 3)
DimPlot(embryos, group.by = "development", reduction = "pca", pt.size = 3)

#Running the pseudotime analysis on all embryos with the genes previously selected
sce_hvg <- subseteR(embryos, dev_genes)
sce_hvg <- reductoR(sce_hvg, embryos, "group")

sce_hvg <- slingshot(sce_hvg, reducedDim = 'PCA')

pt_hvg <- ploteR(sce_hvg, "oc_hvg", embryos, "group")

pt_hvg$stage = embryos$development[rownames(pt_hvg)]
pt_hvg$stage[pt_hvg$stage%in%c("early", "middle", "late", "oocyte")] = "stopped"
pt_hvg$cell = embryos$cell[rownames(pt_hvg)]
pt_hvg$ovulation = embryos$ovulation[rownames(pt_hvg)]

comp_list <- list(c("YS_bad","YS_good")) # bad and good groups correspond to S and SN respectively
res_hvg <- rbind(dif_testeR(pt_hvg[pt_hvg$stage%in%c("morula") & pt_hvg$ovulation == "so",], comp_list),
                 dif_testeR(pt_hvg[pt_hvg$stage%in%c("blastocyst") & pt_hvg$ovulation == "so",], comp_list))

#plot used in Figure 4d
pdf("pseudotime_embryos_withref_withny.pdf", width = 5, height = 4)
par(mar=c(5.1, 8.1, 4.1, 2.1))
maxpt = max(pt_hvg$pseudotime_1, na.rm = T)
maxmoru = max(pt_hvg$pseudotime_1[pt_hvg$stage=="morula"], na.rm = T)
lenwisk = maxpt*0.02
fact_group = factor(pt_hvg$group_label, levels = c("NA", "stopped", "YS_bad", "YS_good", "nat"), labels = c("unclassified", "stopped\ndeveloping", "S", "SN", "NY"))
stripchart(pt_hvg$pseudotime_1 ~ fact_group,
           pch = 20, method = "jitter", xlab = "pseudotime", las = 1, xlim = c(0, maxpt*1.5))
segments(x0 = aggregate(pseudotime_1 ~ group_label + stage, pt_hvg[!pt_hvg$stage%in%c("unclear", "stopped"),], median)$pseudotime_1, 
         y0 = c(1,5,3,4,1,5,3,4)-0.2, 
         y1 = c(1,5,3,4,1,5,3,4)+0.2, col = "red")
segments(x0 = c(maxmoru + lenwisk*3, maxmoru + lenwisk*3, maxmoru + lenwisk*4), y0=c(3,4,3), x1 = rep(maxmoru + lenwisk*4, 3), y1 = c(3,4,4))
segments(x0 = c(maxpt + lenwisk*3, maxpt + lenwisk*3, maxpt + lenwisk*4), y0=c(3,4,3), x1 = rep(maxpt + lenwisk*4, 3), y1 = c(3,4,4))
text(x = maxmoru + lenwisk*4, y = 3.5, label = paste0("pval\n", round(res_hvg$p_value[1], 3)), pos = "4")
text(x = maxpt + lenwisk*4, y = 3.5, label = paste0("pval\n", round(res_hvg$p_value[2], 3)), pos = "4")
dev.off()


#computing correlation as before for comparison
colData_ivf = data.frame(row.names = names(embryos$orig.ident),
                         stage = embryos$group,
                         score = pt_hvg$pseudotime_1[match(names(embryos$orig.ident), rownames(pt_hvg))])

expr_table_ivf = embryos@assays$SCT@counts[rownames(embryos@assays$SCT@counts)%in%dev_genes,]
expr_table_ivf = expr_table_ivf[rowMeans(expr_table_ivf)>1,]

mlm_table_ivf = cbind(colData_ivf, t(expr_table_ivf[,]))
mlm_table_ivf_good = mlm_table_ivf[mlm_table_ivf$stage == "YS_good",]
mlm_table_ivf_bad = mlm_table_ivf[mlm_table_ivf$stage == "YS_bad",]

cor_table_ivf = as.data.frame(t(apply(mlm_table_ivf[,3:ncol(mlm_table_ivf)], 2, function(x){unlist(cor.test(x =x, y = mlm_table_ivf$score, method = "spearman")[c("estimate", "p.value")])})))
cor_table_ivf$gene_name = ensemnl_list_mouse$external_gene_name[match(rownames(cor_table_ivf), ensemnl_list_mouse$ensembl_gene_id)]
cor_table_ivf$padj = p.adjust(cor_table_ivf$p.value, method = "BH")

#All developmental genes have consistent changes in reference and experiment 
#pseudotime analyses (positively or negatively correlated in both)
merge_table = merge(cor_table_val, cor_table_ivf, by = 0)
table(merge_table$estimate.rho.x>0, merge_table$estimate.rho.y>0)
sum(merge_table$estimate.rho.x * merge_table$estimate.rho.y < 0 & merge_table$padj.y<0.05)/nrow(merge_table)*100




