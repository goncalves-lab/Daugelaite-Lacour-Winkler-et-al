#calculating DEGs using pseudobulk approach
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(RColorBrewer)

file_readeR <- function(file_dir){
  file_names <- list.files(file_dir)
  reads_log <- read.table(file_names[1], sep = "\t",fill = T)
  reads_log <- reads_log[5:nrow(reads_log),1:2]
  colnames(reads_log) <- c("gene_names", gsub("_ReadsPerGene.out.tab","",file_names[1]))
  for (i in seq(from=2, to=length(file_names))){
    reads_log_temp <- read.table(file_names[i], sep = "\t",fill = T)
    reads_log_temp <- as.data.frame(reads_log_temp[5:nrow(reads_log_temp),2])
    colnames(reads_log_temp) <- c(gsub("_ReadsPerGene.out.tab","",file_names[i]))
    reads_log <- cbind (reads_log, reads_log_temp)
  }
  return(reads_log)
}

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

thresholdeR <- function(x) {
  if (x[1]<0 & x[2] < 0.05) {
    return("significant_down")
  } 
  if (x[1]>0 & x[2] < 0.05) {
    return("significant_up")
  } 
  else {
    return("non-significant")
  }
}

volcano_plotteR <- function(dge_df, title){
  colnames(dge_df) <- c("gene_name", "mean","log2fc", "padj")
  dge_df=dge_df[!is.na(dge_df$padj),]
  dge_df$label=apply(dge_df[,c('log2fc','padj')], 1, thresholdeR )
  #dge_df$padj=-log10(dge_df$padj)
  #results=results[rev(order(results$label)),]
  # creating color palette
  color1=brewer.pal(9,"Reds")
  #color2=brewer.pal(9,"Purples")
  #color3=brewer.pal(9, "BuGn")
  
  cols <- c("non-significant" = "lightgrey", "significant_up" = "#B31B21", "significant_down"="#1465AC")
  
  # Make a basic ggplot2 object
  vol <- ggplot(dge_df, aes(x = log2(mean), y = log2fc, color = label, labels=gene_name))
  
  # inserting mnaual colors as per color pallette and more
  print(vol +   
          scale_color_manual(values = cols)+
          #ggtitle(title) +
          geom_point(size = 0.4, alpha = 1, na.rm = T)+
          theme_bw(base_size = 14) + 
          theme(legend.position = "right")+
          #geom_hline(yintercept = 1.301029996, colour="#990000", linetype="dashed") +
          #scale_y_continuous(trans = "padj")+
          theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x=expression("log"[2]*" mean"),y=expression("log"[2]*" FC"))+
          theme(axis.text=element_text(size=34), axis.title=element_text(size=36))+
          theme(legend.position="none")+
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+
          theme(plot.title = element_text(hjust = 0.5, size=36))+ ylim(-5,5))
}

volcano_plotteR_inverse <- function(dge_df, title){
  colnames(dge_df) <- c("gene_name", "mean","log2fc", "padj")
  dge_df=dge_df[!is.na(dge_df$padj),]
  dge_df$label=apply(dge_df[,c('log2fc','padj')], 1, thresholdeR )
  #dge_df$padj=-log10(dge_df$padj)
  #results=results[rev(order(results$label)),]
  # creating color palette
  color1=brewer.pal(9,"Reds")
  #color2=brewer.pal(9,"Purples")
  #color3=brewer.pal(9, "BuGn")
  
  cols <- c("non-significant" = "lightgrey", "significant_up" = "#B31B21", "significant_down"="#1465AC")
  dge_df$mean <- -log2(dge_df$mean)
  # Make a basic ggplot2 object
  vol <- ggplot(dge_df, aes(x = mean, y = log2fc, color = label, labels=gene_name))
  
  # inserting mnaual colors as per color pallette and more
  print(vol +   
          scale_color_manual(values = cols)+
          #ggtitle(title) +
          geom_point(size = 0.4, alpha = 1, na.rm = T)+
          theme_bw(base_size = 14) + 
          theme(legend.position = "right")+
          #geom_hline(yintercept = 1.301029996, colour="#990000", linetype="dashed") +
          #scale_y_continuous(trans = "padj")+
          theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x=expression("log"[2]*" mean"),y=expression("log"[2]*" FC"))+
          theme(axis.text=element_text(size=34), axis.title=element_text(size=36))+
          theme(legend.position="none")+
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+
          theme(plot.title = element_text(hjust = 0.5, size=36))+ scale_x_continuous(labels = abs)+
          scale_y_continuous(position = "right", limits = c(-5,5)))
}

#############################################################################

#reading in count data
reads_df <- read.table("reads_all_samples.tsv", sep=",", row.names = 1, header=T)
#reading in metadata
sample_list <- read.csv("metadata_filtered.csv", header = T)
sample_list$Mouse <- as.factor(sample_list$Mouse)
sample_list$ovulation <- as.factor(sample_list$ovulation)
sample_list$age <- factor(sample_list$age, levels = c("Y", "O")) 
reads_df <- reads_df[,colnames(reads_df)%in%sample_list$label]
sample_list <- sample_list[sample_list$label%in%colnames(reads_df),]
######################### 
#performing dge on OCs
########################
sample_list_oc <- sample_list[grep("OC", sample_list$cell_type),]

#summing up reads
summed_df_oc <- data.frame(row.names = row.names(reads_df))
summed_sample_list <- data.frame()
for (i in levels(sample_list$Mouse)){
  samples_df_oc <- sample_list_oc[grep(i, sample_list_oc$Mouse),]
  summed_sample_list <- rbind(summed_sample_list ,samples_df_oc[1,])
  reads_subset <- reads_df[,colnames(reads_df)%in%samples_df_oc$label]
  reads_subset <- as.data.frame(rowSums(reads_subset))
  colnames(reads_subset) <- i
  summed_df_oc <- cbind(summed_df_oc,reads_subset)
}
rownames(summed_sample_list) <- summed_sample_list$Mouse 
summed_sample_list <- summed_sample_list[,2:4]
summed_sample_list$group <- paste0(summed_sample_list$ovulation,"_", summed_sample_list$age)

#condition_list=list(c("NO_Y","NO_O"), c("SO_Y","SO_O"), c("NO_Y","SO_Y"), c("NO_O", "SO_O"))
dge_oc_ov <- analyseR_ovulation(summed_df_oc, summed_sample_list) 
dge_oc_ov$Gene_id <- rownames(dge_oc_ov)
dge_oc_ov <- ensembleR(dge_oc_ov, dge_oc_ov$Gene_id, "Gene_id")
dge_oc_age <- analyseR_age(summed_df_oc, summed_sample_list)
dge_oc_age$Gene_id <- rownames(dge_oc_age)
dge_oc_age <- ensembleR(dge_oc_age, dge_oc_age$Gene_id, "Gene_id")

rld_oc <- normalizeR(summed_df_oc, summed_sample_list)


######################### 
#performing dge on GCs
########################

sample_list_gc <- sample_list[grep("GC", sample_list$cell_type),]
summed_df_gc <- data.frame(row.names = row.names(reads_df))
for (i in levels(sample_list$Mouse)){
  samples_df_gc <- sample_list_gc[grep(i, sample_list_gc$Mouse),]
  reads_subset <- reads_df[,colnames(reads_df)%in%samples_df_gc$label]
  reads_subset <- as.data.frame(rowSums(reads_subset))
  colnames(reads_subset) <- i
  summed_df_gc <- cbind(summed_df_gc,reads_subset)
}


#condition_list=list(c("NO_Y","NO_O"), c("SO_Y","SO_O"), c("NO_Y","SO_Y"), c("NO_O", "SO_O"))
dge_gc_ov <- analyseR_ovulation(summed_df_gc, summed_sample_list) 
dge_gc_ov$Gene_id <- rownames(dge_gc_ov)
dge_gc_ov <- ensembleR(dge_gc_ov, dge_gc_ov$Gene_id, "Gene_id")
dge_gc_age <- analyseR_age(summed_df_gc, summed_sample_list)
dge_gc_age$Gene_id <- rownames(dge_gc_age)
dge_gc_age <- ensembleR(dge_gc_age, dge_gc_age$Gene_id, "Gene_id")

rld_gc <- normalizeR(summed_df_gc, summed_sample_list)


write.csv(dge_oc_age, "dge_oc_age.csv", row.names = F, quote = F)
write.csv(dge_oc_ov, "dge_oc_ov.csv", row.names = F, quote = F)
write.csv(dge_gc_age, "dge_gc_age.csv", row.names = F, quote = F)
write.csv(dge_gc_ov, "dge_gc_ov.csv", row.names = F, quote = F)




