#Differential gene expression analysis using a pseudobulk approach

library(DESeq2)
library(Seurat)
library(biomaRt)


normalizeR <- function(reads, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=reads, colData=metadata, design= ~ group)
  # filtering out genes with 0 counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  # differential expression analysis
  dds <- DESeq(dds)
  rld <- rlog(dds, blind=FALSE)
  return(assay(rld))
}


analyseR_ovulation <- function(reads, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=reads, colData=metadata, design= ~ age+ ovulation + age:ovulation)
  # filtering out genes with 0 counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  # differential expression analysis
  dds <- DESeq(dds)
  res_1 <- as.data.frame(results(dds, contrast=c("ovulation","SO","NO")))
  colnames(res_1) <- c("baseMean_y", "log2FoldChange_y", "lfcSE_y", "stat_y", "pvalue_y","padj_y")
  res_2<- results(dds, list( c("ovulation_SO_vs_NO","ageO.ovulationSO") ))
  colnames(res_2) <- c("baseMean_o", "log2FoldChange_o", "lfcSE_o", "stat_o", "pvalue_o","padj_o")
  res_3 <-results(dds, name="ageO.ovulationSO")
  colnames(res_3) <- c("baseMean_ov", "log2FoldChange_ov", "lfcSE_ov", "stat_ov", "pvalue_ov","padj_ov")
  res <- cbind(res_1, res_2, res_3)
  #merging all result tables in one
  return(res)
}


analyseR_age <- function(reads, metadata) {
  dds <- DESeqDataSetFromMatrix(countData=reads, colData=metadata, design= ~ ovulation + age + ovulation:age)
  # filtering out genes with 0 counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  # differential expression analysis
  dds <- DESeq(dds)
  res_1 <- as.data.frame(results(dds, contrast=c("age","O","Y")))
  colnames(res_1) <- c("baseMean_NO", "log2FoldChange_NO", "lfcSE_NO", "stat_NO", "pvalue_NO","padj_NO")
  res_2<- results(dds, list( c("age_O_vs_Y","ovulationSO.ageO") ))
  colnames(res_2) <- c("baseMean_SO", "log2FoldChange_SO", "lfcSE_SO", "stat_SO", "pvalue_SO","padj_SO")
  res_3 <-results(dds, name="ovulationSO.ageO")
  colnames(res_3) <- c("baseMean_age", "log2FoldChange_age", "lfcSE_age", "stat_age", "pvalue_age","padj_age")
  res <- cbind(res_1, res_2, res_3)
  #merging all result tables in one
  return(res)
}


ensembleR=function(df, df_column, name_column){
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  filters=listFilters(ensembl)
  ensemnl_list=getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values = df_column, mart=ensembl)
  full_list=merge(df, ensemnl_list,by.x=name_column,by.y="ensembl_gene_id")
  return(full_list)
}


#############################################################################

#reading in count data
load("ovulation.Rdata") #Seurat object created by script create_seurat_age_ov.R
reads_df = ovulation@assays$SCT@counts

#reading in metadata
sample_list <- read.csv("metadata_filtered.csv", header = T) #sample metadata table from ArrayExpress
sample_list$Mouse <- as.factor(sample_list$Individual)
sample_list$ovulation <- factor(sample_list$Stimulus, 
                                   levels = c("natural ovulation", "superovulation"),
                                   labels = c("NO", "SO"))
sample_list$age <- factor(sample_list$Age, levels = c(3, 12), labels = c("Y", "O")) 
sample_list$Cell_type <- factor(sample_list$Cell_type, 
                                levels = c("granulosa cell", "oocyte"),
                                labels = c("GC", "OC"))


######################### 
#performing dge on OCs
########################
sample_list_oc <- sample_list[grep("OC", sample_list$Cell_type),]

#summing up reads for pseudobulk
summed_df_oc <- data.frame(row.names = row.names(reads_df))
summed_sample_list <- data.frame()
for (i in levels(sample_list$Mouse)){
  samples_df_oc <- sample_list_oc[grep(i, sample_list_oc$Mouse),] #getting list of samples from each mouse
  summed_sample_list <- rbind(summed_sample_list ,samples_df_oc[1,]) #adding the mouse info to the new metadata
  reads_subset <- reads_df[,colnames(reads_df)%in%samples_df_oc$Name] #subsetting and pseudobulking reads
  reads_subset <- as.data.frame(rowSums(reads_subset))
  colnames(reads_subset) <- i
  summed_df_oc <- cbind(summed_df_oc,reads_subset)
}
rownames(summed_sample_list) <- summed_sample_list$Mouse 
summed_sample_list <- summed_sample_list[,c("Mouse", "ovulation", "age")]
summed_sample_list$group <- paste0(summed_sample_list$ovulation,"_", summed_sample_list$age)

dge_oc_ov <- analyseR_ovulation(summed_df_oc, summed_sample_list) 
dge_oc_ov$Gene_id <- rownames(dge_oc_ov)
dge_oc_ov <- ensembleR(dge_oc_ov, dge_oc_ov$Gene_id, "Gene_id")
dge_oc_age <- analyseR_age(summed_df_oc, summed_sample_list)
dge_oc_age$Gene_id <- rownames(dge_oc_age)
dge_oc_age <- ensembleR(dge_oc_age, dge_oc_age$Gene_id, "Gene_id")

#rld_oc <- normalizeR(summed_df_oc, summed_sample_list)

######################### 
#performing dge on GCs
########################

#summing up reads for pseudobulk
sample_list_gc <- sample_list[grep("GC", sample_list$Cell_type),]
summed_df_gc <- data.frame(row.names = row.names(reads_df))
for (i in levels(sample_list$Mouse)){
  samples_df_gc <- sample_list_gc[grep(i, sample_list_gc$Mouse),]
  reads_subset <- reads_df[,colnames(reads_df)%in%samples_df_gc$Name]
  reads_subset <- as.data.frame(rowSums(reads_subset))
  colnames(reads_subset) <- i
  summed_df_gc <- cbind(summed_df_gc,reads_subset)
}


dge_gc_ov <- analyseR_ovulation(summed_df_gc, summed_sample_list) 
dge_gc_ov$Gene_id <- rownames(dge_gc_ov)
dge_gc_ov <- ensembleR(dge_gc_ov, dge_gc_ov$Gene_id, "Gene_id")
dge_gc_age <- analyseR_age(summed_df_gc, summed_sample_list)
dge_gc_age$Gene_id <- rownames(dge_gc_age)
dge_gc_age <- ensembleR(dge_gc_age, dge_gc_age$Gene_id, "Gene_id")

#rld_gc <- normalizeR(summed_df_gc, summed_sample_list)


write.csv(dge_oc_age, "dge_oc_age.csv", row.names = F, quote = F)
write.csv(dge_oc_ov, "dge_oc_ov.csv", row.names = F, quote = F)
write.csv(dge_gc_age, "dge_gc_age.csv", row.names = F, quote = F)
write.csv(dge_gc_ov, "dge_gc_ov.csv", row.names = F, quote = F)


