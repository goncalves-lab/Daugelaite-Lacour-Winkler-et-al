#performs ora analysis
library(data.table)
library(stringr)
library(biomaRt)
library(RColorBrewer)
library(tidyverse)

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


thresholdeR <- function(x) {
  if (abs(x[1]) == Inf && !is.na(x[1])) {
    return(0)
  } 
  else {
    return(x[2])
  }
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
    num_vector_1=c((length(in_set_de_v)-1),length(gene_ds$gene_name),length(gene_nds$gene_name), (length(in_set_de_v)+length(in_set_nde_v)))
    #num_vector_2=c(length(in_set_de_v),length(df_1_de$entrezgene_id)-length(in_set_de_v),length(in_set_nde_v), length(df_1_nde$entrezgene_id)-length(in_set_nde_v))
    p_value_1=phyper(num_vector_1[1], num_vector_1[2], 21859-num_vector_1[2], num_vector_1[4], lower.tail = FALSE, log.p = FALSE)
    #p_value_2=fisher.test(matrix(c(num_vector_2[1], num_vector_2[2], num_vector_2[3], num_vector_2[4]), 2, 2), alternative='less')$p.value
    num_vector_1[5]=p_value_1
    df_hyper=rbind(df_hyper,c(gs_mmu_data[[i]][[1]],gs_mmu_data[[i]][[2]],
                              unname(as.list(num_vector_1))))
    df_hyper[,1]=as.character(df_hyper[,1])
    df_hyper[,2]=as.character(df_hyper[,2])
  }
  colnames(df_hyper)=c("pathway","description","in_set_df","df","n_df","in_set","p_value")
  df_hyper$p_value_adj=p.adjust(df_hyper$p_value, method = "BH", n = length(df_hyper$p_value))
  df_hyper$in_set_df=df_hyper$in_set_df+1
  df_hyper$gene_ratio=df_hyper$in_set_df/df_hyper$in_set
  write.table(df_hyper,file = file_name, quote = F,row.names = F,sep=",")
  return(df_hyper)
}

enricheR_common=function(total_genes, sig_genes, gs_mmu_data, file_name){
  df_hyper <- data.frame()
  for (i in seq(from=1, to =length(gs_mmu_data))){
    in_set_de_v=sig_genes%in%gs_mmu_data[[i]][[3]]
    in_set_de_v=in_set_de_v[in_set_de_v==TRUE]
    num_vector_1=c((length(in_set_de_v)-1),length(sig_genes),total_genes-length(sig_genes), length(gs_mmu_data[[i]][[3]]))
    #num_vector_2=c(length(in_set_de_v),length(df_1_de$entrezgene_id)-length(in_set_de_v),length(in_set_nde_v), length(df_1_nde$entrezgene_id)-length(in_set_nde_v))
    p_value_1=phyper(num_vector_1[1], num_vector_1[2], 21859-num_vector_1[2], num_vector_1[4], lower.tail = FALSE, log.p = FALSE)
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
  write.table(df_hyper,file = file_name, quote = F,row.names = F,sep=",")
  return(df_hyper)
}



col_1 = colorRampPalette(rev(brewer.pal(9, "RdBu")))
plotteR=function(df_gse, subset_val){
  df_gse <- df_gse[order(df_gse$p_value_adj),]
  print(ggplot(df_gse[subset_val,], aes(gene_ration, pathway))+
          geom_point(aes(colour=p_value_adj, size=in_set_df))+
          scale_colour_gradientn(colours = col_1(200), name="padj value")+
          theme_classic()+
          theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
          labs(x="Gene ratio", y="",size="")+
          theme(legend.position="right", 
                legend.title = element_text(colour="Black", size=18),
                legend.text = element_text(colour="black", size = 20))+
          scale_size(range = c(1,25)))
}

ensembleR=function(val,attr,filt){
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  filters=listFilters(ensembl)
  attributes=listAttributes(ensembl)
  ensembl_df=getBM(attributes=attr,filters=filt,values = val,mart=ensembl)
  return(ensembl_df)
}

h_mmu_data=GeneSetListFromGmt("h_all_v7.1_symbols_mm.gmt")
c2_mmu_data=GeneSetListFromGmt("c2_cp_v7.1_symbols_mm.gmt")
c5bp_mmu_data=GeneSetListFromGmt("c5_bp_v7.1_symbols_mm.gmt")

#produced in script ora.R
dge_oc_ov <- read.csv("dge_oc_ov.csv", header=T)
dge_oc_age <- read.csv("dge_oc_age.csv", header=T)

dge_oc_ov_nsig <- dge_oc_ov[dge_oc_ov$padj_y>0.05,]
dge_oc_age_nsig <- dge_oc_age[dge_oc_age$padj_NO>0.05,]
dge_oc_ov_nsig <- dge_oc_ov_nsig[!is.na(dge_oc_ov_nsig$Gene_id),]
dge_oc_age_nsig <- dge_oc_age_nsig[!is.na(dge_oc_age_nsig$Gene_id),]
dge_oc_nsig <- intersect(dge_oc_age_nsig$external_gene_name, dge_oc_ov_nsig$external_gene_name)

dge_oc_ov_sig <- dge_oc_ov[dge_oc_ov$padj_y<0.05,]
dge_oc_age_sig <- dge_oc_age[dge_oc_age$padj_NO<0.05,]
dge_oc_ov_sig <- dge_oc_ov_sig[!is.na(dge_oc_ov_sig$Gene_id),]
dge_oc_age_sig <- dge_oc_age_sig[!is.na(dge_oc_age_sig$Gene_id),]
dge_oc_sig <- intersect(dge_oc_age_sig$external_gene_name, dge_oc_ov_sig$external_gene_name)

dge_oc_sig <- data.frame(p_value=0.01, gene=dge_oc_sig)
dge_oc_nsig <- data.frame(p_value=1, gene=dge_oc_nsig)
dge_oc <- rbind(dge_oc_sig, dge_oc_nsig)


h_y_ov <- enricheR(dge_oc, "h_so_age_oc_merged.csv", h_mmu_data)
c2_y_ov <- enricheR(dge_oc, "c2_so_age_oc_merged.csv", c2_mmu_data)
c5bp_y_ov <-enricheR(dge_oc,  "c5bp_so_age_oc_merged.csv", c5bp_mmu_data)

y_ov_oc <- bind_rows(h_y_ov, c2_y_ov, c5bp_y_ov)
y_ov_oc$p_value_adj_final <- p.adjust(y_ov_oc$p_value, method = "BH")

dge_gc_ov <- read.csv("dge_gc_ov.csv", header=T)
dge_gc_age <- read.csv("dge_gc_age.csv", header=T)

dge_gc_ov_nsig <- dge_gc_ov[dge_gc_ov$padj_y>0.05,]
dge_gc_age_nsig <- dge_gc_age[dge_gc_age$padj_NO>0.05,]
dge_gc_ov_nsig <- dge_gc_ov_nsig[!is.na(dge_gc_ov_nsig$Gene_id),]
dge_gc_age_nsig <- dge_gc_age_nsig[!is.na(dge_gc_age_nsig$Gene_id),]
dge_gc_nsig <- intersect(dge_gc_age_nsig$external_gene_name, dge_gc_ov_nsig$external_gene_name)

dge_gc_ov_sig <- dge_gc_ov[dge_gc_ov$padj_y<0.05,]
dge_gc_age_sig <- dge_gc_age[dge_gc_age$padj_NO<0.05,]
dge_gc_ov_sig <- dge_gc_ov_sig[!is.na(dge_gc_ov_sig$Gene_id),]
dge_gc_age_sig <- dge_gc_age_sig[!is.na(dge_gc_age_sig$Gene_id),]
dge_gc_sig <- intersect(dge_gc_age_sig$external_gene_name, dge_gc_ov_sig$external_gene_name)

dge_gc_sig <- data.frame(p_value=0.01, gene=dge_gc_sig)
dge_gc_nsig <- data.frame(p_value=1, gene=dge_gc_nsig)
dge_gc <- rbind(dge_gc_sig, dge_gc_nsig)


h_y_ov_gc <- enricheR(dge_gc, "h_so_age_gc_merged.csv", h_mmu_data)
c2_y_ov_gc <- enricheR(dge_gc, "c2_so_age_gc_merged.csv", c2_mmu_data)
c5bp_y_ov_gc <-enricheR(dge_gc,  "c5bp_so_age_gc_merged.csv", c5bp_mmu_data)

y_ov_gc <- bind_rows(h_y_ov_gc, c2_y_ov_gc, c5bp_y_ov_gc)
y_ov_gc$p_value_adj_final <- p.adjust(y_ov_gc$p_value, method = "BH")

toMatch_oc <- toupper(c("meiosis", "spindle", "checkpoint",  "Chromosome segregation",
"Estrogen", "insulin", "PI3Akt", "MAPK", "GONADOTROPIN", "DNA repair", "oxidative", "apoptosis", "STEROID", "PROGESTERONE", "ANEUPLOIDY"))
oc_pathways <- y_ov_oc[y_ov_oc$pathway%in%(grep(paste(toMatch_oc,collapse="|"), y_ov_oc$pathway, value=TRUE)),]
oc_pathways <- oc_pathways[oc_pathways$p_value_adj_final<0.05,]
write.csv(oc_pathways, "oc_pathways_common.csv", row.names = F, quote = F)

toMatch_gc <- toupper(c("Proliferation", "cell_cycle", "Mitosis","G2M", "Steroidogenesis", "Estrogen", "steroid", "extracellular",
                        "repair", "Apoptosis", "ANDROGEN", "PROGESTERONE",
                        "FATTY", "LUTEIN", "MAPK", "KIT","ERK", "EGF", "EPIDERMAL GROWTH"))
gc_pathways <- y_ov_gc[y_ov_gc$pathway%in%(grep(paste(toMatch_gc,collapse="|"), y_ov_gc$pathway, value=TRUE)),]
gc_pathways <- gc_pathways[gc_pathways$p_value_adj_final<0.05,]
write.csv(gc_pathways, "gc_pathways_common.csv", row.names = F, quote = F)
