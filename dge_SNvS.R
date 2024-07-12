#Selecting genes for classifier

#-----Library--------
library(Seurat)
library(biomaRt)
library(RColorBrewer)
library(fgsea)
library(stringr)
library(DESeq2)
library(ggplot2)

#-----Functions--------

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
enricheR=function(gene_ds_total, file_name, gs_mmu_data){
  df_hyper <- data.frame()
  colnames(gene_ds_total) <- c("padj", "gene_name")
  gene_ds <- gene_ds_total[gene_ds_total$padj<0.05,]
  gene_ds <- gene_ds[!is.na(gene_ds$padj),]
  gene_nds <- gene_ds_total[gene_ds_total$padj>0.05,]
  gene_nds <- gene_nds[!is.na(gene_nds$padj),]
  for (i in seq(from=1, to =length(gs_mmu_data))){
    in_set_de_v=as.character(gene_ds$gene_name)%in%gs_mmu_data[[i]][[3]]
    in_set_de_v=in_set_de_v[in_set_de_v==TRUE]
    in_set_nde_v=as.character(gene_nds$gene_name)%in%gs_mmu_data[[i]][[3]]
    in_set_nde_v=in_set_nde_v[in_set_nde_v==TRUE]  
    num_vector_1=c((length(in_set_de_v)-1),length(gene_ds$gene_name),length(gene_nds$gene_name), length(in_set_de_v)+length(in_set_nde_v))
    #num_vector_2=c(length(in_set_de_v),length(df_1_de$entrezgene_id)-length(in_set_de_v),length(in_set_nde_v), length(df_1_nde$entrezgene_id)-length(in_set_nde_v))
    p_value_1=phyper(num_vector_1[1], num_vector_1[2], num_vector_1[3], num_vector_1[4], lower.tail = FALSE, log.p = FALSE)
    #p_value_2=fisher.test(matrix(c(num_vector_2[1], num_vector_2[2], num_vector_2[3], num_vector_2[4]), 2, 2), alternative='less')$p.value
    num_vector_1[5]=p_value_1
    df_hyper=rbind(df_hyper,c(gs_mmu_data[[i]][[1]],gs_mmu_data[[i]][[2]],unname(as.list(num_vector_1))))
    df_hyper[,1]=as.character(df_hyper[,1])
    df_hyper[,2]=as.character(df_hyper[,2])
  }
  colnames(df_hyper)=c("pathway","description","in_set_df","df","n_df","in_set","p_value")
  df_hyper$p_value_adj=p.adjust(df_hyper$p_value, method = "BH", n = length(df_hyper$p_value))
  df_hyper$in_set_df=df_hyper$in_set_df+1
  df_hyper$gene_ration=df_hyper$in_set_df/df_hyper$df
  #write.table(df_hyper,file = file_name, quote = F,row.names = F,sep=",")
  return(df_hyper)
}


permutation_test = function(X, label) {
  t_dist = data.frame()
  level1 = levels(label)[1]
  level2 = levels(label)[2]
  real_FC = rowMeans(X[,label==level1]) / rowMeans(X[,label==level2])
  for (i in 1:1000) {
    new_label = sample(label, length(label), replace = F)
    t_dist = rbind(t_dist, rowMeans(X[,new_label==level1]) / rowMeans(X[,new_label==level2]))
  }
  pval_df = data.frame(row.names = rownames(X), pval = rep(NA, nrow(X)))
  for (i in 1:ncol(t_dist)) {
    gene = rownames(X)[i]
    pval_df[gene, "pval"] = ecdf(t_dist[,i])(real_FC[gene])
  }
  pval_df$padj = p.adjust(pval_df$pval)
  return(pval_df)
}


#-----Parameters and data--------

pathways_h <- gmtPathways("gene_sets/h.all.v7.2.symbols.gmt")
pathways_c2 <- gmtPathways("gene_sets/c2.cp.reactome.v7.2.symbols.gmt")
pathways_go = gmtPathways("gene_sets/c5.go.bp.v7.2.symbols.gmt")

ensembl_mouse=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list_mouse=getBM(attributes=c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_associated_gene_name", "chromosome_name"), mart=ensembl_mouse)


load("clusters_Y_allLR.Rdata") # clusters obtained from cell communication analysis


load("ovulation.Rdata") #Seurat object

#matching OC and GC labels
paired_id = data.frame()
gc_list = names(ovulation$label)[grep("GC",names(ovulation$label))]
oc_list = names(ovulation$label)[grep("OC",names(ovulation$label))]
for (i in 1:length(gc_list)) {
  if (length(grep(gsub("GC","OC",gc_list[i]),oc_list))) {
    paired_id = rbind(paired_id, c(gc_list[i], oc_list[grep(gsub("GC","OC",gc_list[i]),oc_list)]))
  }
}
colnames(paired_id) = c("GC", "OC")
paired_id$group = paste(ovulation$age[paired_id$GC],ovulation$ovulation[paired_id$GC],sep="_")

#adding information about clustering in cell comm issues
paired_id$group_cellcom = paired_id$group
paired_id$group_cellcom[paired_id$OC%in%names(clusters_cell_comm)[clusters_cell_comm=="YS_good"]] = "Y_SO_good"
paired_id$group_cellcom[paired_id$OC%in%names(clusters_cell_comm)[clusters_cell_comm=="YS_bad"]] = "Y_SO_bad"



#-------Main--------

dds <- DESeqDataSetFromMatrix(countData = ovulation@assays$RNA[,paired_id$GC[paired_id$group%in%c("Y_SO")]],
                              colData = data.frame(row.names = paired_id$GC[paired_id$group%in%c("Y_SO")],
                                                   group = paired_id$group_cellcom[paired_id$group%in%c("Y_SO")],
                                                   mouse = ovulation$mouse[paired_id$GC[paired_id$group%in%c("Y_SO")]]),
                              design= ~ group + mouse)
dds_groupmouse <- DESeq(dds)
normcounts = counts(dds_groupmouse, normalized=TRUE)
res_dds_groupmouse <- as.data.frame(results(dds_groupmouse, name="group_Y_SO_good_vs_Y_SO_bad"))
res_dds_groupmouse$genename=ensemnl_list_mouse$external_gene_name[match(rownames(res_dds_groupmouse),ensemnl_list_mouse$ensembl_gene_id)]


# res_dds_groupmouse$humangene = ensemnl_list_mouse$hsapiens_homolog_associated_gene_name[match(rownames(res_dds_groupmouse), ensemnl_list_mouse$ensembl_gene_id)]
# pathways_enriched = enricheR(res_dds_groupmouse[,c("padj", "humangene")], 
#                              gs_mmu_data = GeneSetListFromGmt("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/gene_sets/h.all.v7.2.symbols.gmt"))
# pathways_enriched = enricheR(res_dds_groupmouse[,c("padj", "humangene")], 
#                              gs_mmu_data = GeneSetListFromGmt("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/gene_sets/c5.go.bp.v7.2.symbols.gmt"))





gr = c("Y_SO", "Y_NO", "O_NO")
#Selecting genes based on significance, logFoldChange and expression level
signifgenes = rownames(res_dds_groupmouse)[res_dds_groupmouse$padj<0.01 & 
                                             !is.na(res_dds_groupmouse$padj) & 
                                             res_dds_groupmouse$baseMean>100 & 
                                             abs(res_dds_groupmouse$log2FoldChange) > 0.58]
names(signifgenes) = ensemnl_list_mouse$external_gene_name[match(signifgenes, ensemnl_list_mouse$ensembl_gene_id)]
GC_SO_signifgenes = ovulation@assays$SCT@counts[signifgenes, paired_id$GC[paired_id$group%in%gr]]
pca_signifgenes = prcomp(t(GC_SO_signifgenes), scale. = T)
plot(pca_signifgenes$x[,1], pca_signifgenes$x[,2], pch = 20,
     col = as.character(factor(paired_id$group_cellcom[paired_id$group%in%gr], 
                               levels = c("Y_NO",  "Y_SO_bad", "Y_SO_good", "O_NO"), 
                               labels = c("black", "red", "green", "grey"))))
#Since PC1 is separating the two groups, filtering based on link between gene and PC1
signifgenes_pca_filtered = signifgenes[abs(pca_signifgenes$rotation[,1])>0.1]

write.table(res_dds_groupmouse[signifgenes_pca_filtered,], file = "quality_marker_genes.tsv",
            quote = F, sep = ",", row.names = T)

