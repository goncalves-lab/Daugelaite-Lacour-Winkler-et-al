# Differential gene expression analysis for human GC data

#-------Libraries----------

library(DESeq2)
library(biomaRt)
library(fgsea)

#-------Functions---------

normalizeR = function(reads, metadata) {
  dds = DESeqDataSetFromMatrix(countData=reads, colData=metadata, design= ~ group)
  # filtering out genes with 0 counts
  dds = dds[rowSums(counts(dds)) > 1, ]
  # differential expression analysis
  dds = DESeq(dds)
  rld = rlog(dds, blind=FALSE)
  return(assay(rld))
}


ensembleR=function(df,val,attr,filt,x_merge,y_merge){
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  attributes=listAttributes(ensembl)
  ensemnl_list=getBM(attributes=attr,filters=filt,values = val,mart=ensembl)
  full_list=merge(df,ensemnl_list,by.x=x_merge,by.y=y_merge)
  return(full_list)
}

enricheR = function(dge_df,gs_mmu_data){ #preparing data and running fgsea 
  colnames(dge_df) = c("external_gene_name", "stat", "padj", "hsapiens_homolog_associated_gene_name")
  dge_df = dge_df[!duplicated(dge_df$hsapiens_homolog_associated_gene_name),]
  dge_df = dge_df[dge_df$hsapiens_homolog_associated_gene_name!="",]
  dge_df = dge_df[!is.na(dge_df$padj),]
  dge_df = dge_df[order(dge_df$stat),]
  df_t_stat = dge_df$stat
  names(df_t_stat)=dge_df$hsapiens_homolog_associated_gene_name
  fgseaRes = fgsea(gs_mmu_data, df_t_stat, minSize=15, maxSize=500)
  return(fgseaRes)
}    

#-------Data----------

pathways_h = gmtPathways("gene_sets/h.all.v7.2.symbols.gmt")
pathways_c2 = gmtPathways("gene_sets/c2.cp.reactome.v7.2.symbols.gmt")
pathways_go = gmtPathways("gene_sets/c5.go.bp.v7.2.symbols.gmt")

#mouse data
dge_gc = read.csv("dge_gc_age.csv", header = T) #output of the script dge.R
dge_gc = ensembleR(dge_gc, dge_gc$external_gene_name,c("external_gene_name","hsapiens_homolog_associated_gene_name"),"external_gene_name","external_gene_name","external_gene_name")

pathways_all = c(pathways_c2,
                 pathways_h,
                 pathways_go)


#Setting the parameters
predictor = "group"
pathways_to_test = pathways_all

#human data
reads_df=read.table(file = "count_table_all_samples.txt", header = T, row.names = 1) #Available on ArrayExpress

#reading in metadata
sample_list = read.csv("sample_table.csv") #sample metadata table from ArrayExpress
sample_list = sample_list[, c("Name", "Age")]
sample_list$group = factor(sample_list$Age < 35, levels = c(T, F), labels = c("Y", "O"))
sample_list$group[sample_list$Age>32 & sample_list$Age<40] = NA
colnames(sample_list) = c("patient", "age", "group")
sample_list = sample_list[sample_list$patient%in%colnames(reads_df),]
sample_list = sample_list[!duplicated(sample_list$patient),]

#gene_id and names
gene_names = cbind(gene_id = rownames(reads_df), name = reads_df$gene_name)

#-------Main-------------

reads_df = reads_df[, sample_list$patient]

sample_list$predictor = sample_list[, predictor]


dds = DESeqDataSetFromMatrix(countData=reads_df, colData=sample_list, design= ~ predictor)
# filtering out genes with 0 counts
dds = dds[rowSums(counts(dds)) > 1, ]
# differential expression analysis
dds = DESeq(dds)

res = as.data.frame(results(dds))
res$gene_name = gene_names$name[match(rownames(res), gene_names$gene_id)]



## comparison with mouse data

dge_gc_so_y_o = dge_gc[ , c("external_gene_name", "stat_SO", "padj_SO","hsapiens_homolog_associated_gene_name")]
gsea_gc_so_y_o = enricheR(dge_gc_so_y_o, pathways_to_test)


dge_human_age = res[,c("gene_name", "stat", "padj", "gene_name")]
gsea_human = enricheR(dge_human_age, pathways_to_test)

#merging the results from both species
gsea_df = data.frame(pathway = names(pathways_to_test), 
                     human = gsea_human$NES[match(names(pathways_to_test), gsea_human$pathway)],
                     human_padj = gsea_human$padj[match(names(pathways_to_test), gsea_human$pathway)],
                     mouse_SO = gsea_gc_so_y_o$NES[match(names(pathways_to_test), gsea_gc_so_y_o$pathway)],
                     mouse_SO_padj = gsea_gc_so_y_o$padj[match(names(pathways_to_test), gsea_gc_so_y_o$pathway)])

gsea_df = na.exclude(gsea_df[gsea_df$human_padj<0.05 | gsea_df$mouse_SO_padj<0.05,])

write.csv(gsea_df, file = paste0("gsea_human_mouse.csv"))

