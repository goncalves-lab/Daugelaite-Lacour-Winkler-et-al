# Comparing total RNA-seq and SMART-seq2 results

#------Libraries---------

library(Seurat)
library(biomaRt)



#------Data--------------

ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list=getBM(attributes=c("ensembl_gene_id","external_gene_name", "gene_biotype"), mart=ensembl)

#total RNA data
load("tot_RNA.RData") # total RNA-seq Seurat object made by the script create_seurat_totalRNA.R

reads_tot = tot_RNA@assays$SCT@counts
rownames(reads_tot) = sapply(rownames(reads_tot), function(x){strsplit(x, "-")[[1]][2]})


#smart-seq2 data
load("ovulation.Rdata") # SMART-seq2 Seurat object made by the script create_seurat_age_ov.R


#Changing to gene name instead of ensemblID
reads_smart = as.matrix(ovulation@assays$SCT@counts)
reads_smart = reads_smart[, ovulation$cell == "OC" & ovulation$age == "Y"]
rownames(reads_smart) = ensemnl_list$external_gene_name[match(rownames(reads_smart), ensemnl_list$ensembl_gene_id)]


#gene list
polyA = c("Dazl", "Cnot7", "Btg4", "Oosp1", "Oosp3", "Tcl1", "Obox5", "Cpeb4")
degraded = c("Cpeb1", "Nlrp14", "Padi6", "Zp2", "Zp3")
deA = c("Smc4", "Zp1", "Oosp2")

#------Main-----------

log2FC_smart = log2(rowMeans(reads_smart[,grepl("NO", colnames(reads_smart))]) / rowMeans(reads_smart[,grepl("SO", colnames(reads_smart))]))
log2FC_tot = log2(rowMeans(reads_tot[,grepl("NO", colnames(reads_tot))]) / rowMeans(reads_tot[,grepl("SO", colnames(reads_tot))]))

#Selecting a few genes that do not respond to the treatment in both technologies
constant_genes_smart = names(log2FC_smart)[abs(log2FC_smart)<0.005 & is.finite(log2FC_smart)]
constant_genes_tot = names(log2FC_tot)[abs(log2FC_tot)<0.01 & is.finite(log2FC_tot)]
constant_genes = intersect(constant_genes_smart, constant_genes_tot)

#Selecting a few genes that are upregulated by the treatment in both technologies
up_genes_smart = dge_oc_ov_smart$Gene_id[dge_oc_ov_smart$padj_y < 0.05 & !is.na(dge_oc_ov_smart$padj_y) & dge_oc_ov_smart$log2FoldChange_y < 0]
up_genes_tot = dge_oc_ov_tot$Gene_id[dge_oc_ov_tot$padj < 0.1 & !is.na(dge_oc_ov_tot$padj) & dge_oc_ov_tot$log2FoldChange < 0]
up_genes = intersect(up_genes_smart, up_genes_tot)
up_genes = ensemnl_list$external_gene_name[match(up_genes, ensemnl_list$ensembl_gene_id)]



# Extracting the genes of interest in SMART-seq table
matr_smart = reads_smart[c(polyA, deA, degraded, up_genes, constant_genes), ]

# For each cell in one group, we compute the log2 fold change relative to the average of the other group
mean_SO = rowMeans(matr_smart[, grepl("SO", colnames(matr_smart))]) + 1
mean_NO = rowMeans(matr_smart[, grepl("NO", colnames(matr_smart))]) + 1
norm_matr_smart = cbind(log2((1+matr_smart[, grepl("NO", colnames(matr_smart))])/mean_SO),
                        log2((1+matr_smart[, grepl("SO", colnames(matr_smart))])/mean_NO))

# Extracting the genes of interest in total RNA-seq table
matr_tot = reads_tot[c(polyA, deA, degraded, up_genes, constant_genes), ]

# For each cell in one group, we compute the log2 fold change relative to the average of the other group
mean_SO = rowMeans(matr_tot[, grepl("SO", colnames(matr_tot))]) + 1
mean_NO = rowMeans(matr_tot[, grepl("NO", colnames(matr_tot))]) + 1
norm_matr_tot = cbind(log2((1+matr_tot[, grepl("NO", colnames(matr_tot))])/mean_SO),
                      log2((1+matr_tot[, grepl("SO", colnames(matr_tot))])/mean_NO))

norm_matr = cbind(norm_matr_smart, norm_matr_tot)

#Computing graphical parameters
gaps_row = cumsum(c(length(polyA), length(deA), length(degraded), length(up_genes), length(constant_genes)))
gaps_col = ncol(norm_matr_smart)

anno = data.frame(row.names = colnames(norm_matr), ov = substr(colnames(norm_matr), 1,2))
anno_color = list(ov = c("SO"= "#727272", "NO" = "black"))


pheatmap::pheatmap(norm_matr, scale = "none", cluster_cols = F, cluster_rows = F, 
                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu")[11:1])(100), breaks = seq(-2,2,length.out=100),
                   annotation_col = anno, gaps_row = gaps_row, gaps_col = gaps_col,
                   border_color = NA, annotation_colors = anno_color, annotation_names_col = F, 
                   show_colnames = F)


#to subset the heatmap
set.seed(123)
n_smartN = length(grep("NO", colnames(matr_smart)))
n_smartS = length(grep("SO", colnames(matr_smart)))
n_totN = length(grep("NO", colnames(matr_tot)))
n_totS = length(grep("SO", colnames(matr_tot)))
n_groups = cumsum(c(n_smartN, n_smartS, n_totN, n_totS))
samps = c(sample(rownames(anno)[1:n_groups[1]], 10, replace = F),
          sample(rownames(anno)[(n_groups[1]+1):n_groups[2]], 10, replace = F),
          sample(rownames(anno)[(n_groups[2]+1):n_groups[3]], 10, replace = F),
          sample(rownames(anno)[(n_groups[3]+1):n_groups[4]], 10, replace = F))

norm_matr_sub = norm_matr[, samps]
gaps_row = cumsum(c(length(polyA), length(deA), length(degraded), length(up_genes), length(constant_genes)))
gaps_col = 20

pheatmap::pheatmap(norm_matr_sub, scale = "none", cluster_cols = F, cluster_rows = F, 
                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu")[11:1])(100), breaks = seq(-2,2,length.out=100),
                   annotation_col = anno, gaps_row = gaps_row, gaps_col = gaps_col,
                   border_color = NA, annotation_colors = anno_color, annotation_names_col = F, 
                   show_colnames = F)
