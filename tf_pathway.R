library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(ggalluvial)
library(ComplexHeatmap)
library(data.table)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

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

pathway_prep <- function(pathways, pathways_selected){
  pathways_df <- data.frame()
  for (i in seq(from=1, to =length(pathways))){
    if (pathways[[i]][[1]]%in%pathways_selected){
      pathways_name_temp <- pathways[[i]][[1]]
      pathways_gene_temp <- paste(pathways[[i]][[3]], collapse=" ")
      pathways_df_temp <- data.frame(pathway=pathways_name_temp, genes=pathways_gene_temp)
      pathways_df <- rbind(pathways_df_temp, pathways_df)
    }
  }
  return(pathways_df)
}

#hypergeometric test
enricheR <- function(contrast_var, sig_filter, c_pw, regulons, pathways_df){
  c_pw <- c_pw[c_pw$Comparison==contrast_var,]
  if (sig_filter){
    c_pw <- c_pw[c_pw$p_adj<0.05,]}
  regulons$tf <- gsub("(+)", "", regulons$tf, fixed=T)
  c_pw$tf <- substr(c_pw$tf, 1, nchar(c_pw$tf)-3)
  regulons <- regulons[regulons$tf%in%c_pw$tf,]
  bk <- 21833
  tf_l_df <- data.frame(filler=1:nrow(pathways_df))
  tf <- c()
  for(i in seq(from=1, to=nrow(regulons))){
    tf_targets <- strsplit(regulons[i,2], " ")[[1]]
    pathways_gs <- strsplit(pathways_df$genes, " ")
    tf_l <- c()
    tf <- c(regulons[i,1],tf)
    for (j in seq(from=1, to=length(pathways_gs))){
      pathways_gs_s <- pathways_gs[[j]]
      pathways_gs_s <- toupper(pathways_gs_s)
      tf_targets <- toupper(tf_targets)
      test_res <- phyper(length(intersect(tf_targets, pathways_gs_s))-1, length(tf_targets), bk-length(tf_targets) ,length(pathways_gs_s), lower.tail = F, log.p = F)
      tf_l <- c(tf_l, test_res)
    }
    
    names(tf_l) <- pathways_df[,1]
    tf_l <- p.adjust(tf_l, method = "BH")
    tf_l_df <- cbind( tf_l, tf_l_df)
    rownames(tf_l_df) <- names(tf_l)
  }
  
  
  colnames(tf_l_df) <- tf
  tf_l_df <- tf_l_df[,1:ncol(tf_l_df)-1]
  return(tf_l_df)
}

JI_calculatoR <- function(pathways_df, regulons){
  tf_l_ji <- data.frame(filler=1:nrow(pathways_df))
  for(i in seq(from=1, to=nrow(regulons))){
    tf_targets <- strsplit(regulons[i,2], " ")[[1]]
    ligand_gs <- strsplit(pathways_df$genes, " ")
    ji_f <- c()
    for (j in seq(from=1, to=length(ligand_gs))){
      ligand_gs_s <- ligand_gs[[j]]
      ligand_gs_s <- toupper(ligand_gs_s)
      tf_targets <- toupper(tf_targets)
      ji <- length(intersect(tf_targets, ligand_gs_s))/length(union(tf_targets, ligand_gs_s))
      ji_f <- c(ji_f, ji)
    }
    
    names(ji_f) <- pathways_df[,1]
    tf_l_ji <- cbind(tf_l_ji, ji_f)
    rownames(tf_l_ji) <- names(ji_f)
  }
  
  tf_l_ji <- tf_l_ji[,2:ncol(tf_l_ji)]
  colnames(tf_l_ji) <- regulons$tf
  return(tf_l_ji)
  
}

tf_selectoR <- function(tf_l_df, tf_l_ji){
  tf_list <- list()
  for (i in seq(nrow(tf_l_df))){
    print(i)
    tf_p <- t(tf_l_df[i,])
    pathway_name <- colnames(tf_p)
    tf_j <- t(tf_l_ji[i,])
    rownames(tf_j) <- substr(rownames(tf_j,1),1, nchar(rownames(tf_j))-3)
    tf <- merge(tf_p, tf_j, by=0)
    colnames(tf) <- c("tf", "p_adj", "ji")
    tf <-tf[order(tf$p_adj),]
    tf <- tf[tf$p_adj<0.01,]
    tf_list_t <- list(tf$tf)
    names(tf_list_t) <- pathway_name
    tf_list <- c(tf_list_t, tf_list)
  }
  
  tf_list <- Filter(length, tf_list)
  tf_df_all <- data.frame()
  for (i in seq(length(tf_list))){
    print(i)
    tf_df <- data.frame(tf=tf_list[[i]], pathway=names(tf_list)[i])
    tf_df_all <- rbind(tf_df, tf_df_all)
  }
  return(tf_df_all)
}
#####################################################################################################################
#pathway prep
h_mmu_data=GeneSetListFromGmt("h_all_v7.1_symbols_mm.gmt")
c2_mmu_data=GeneSetListFromGmt("c2_cp_v7.1_symbols_mm.gmt")
c5bp_mmu_data=GeneSetListFromGmt("c5_bp_v7.1_symbols_mm.gmt")

pathways <- c(h_mmu_data, c2_mmu_data, c5bp_mmu_data)
regulons_gc <- read.csv("reg_gc") #tfs and their targets
c_pw_gc <-loadRData("gc_separate_pw_nl.Rdata") #produced in script scenic_post.R

pathways_selected <- c("HALLMARK_DNA_REPAIR", "REACTOME_PROTEIN_REPAIR", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                       "GO_POSTREPLICATION_REPAIR", "REACTOME_FATTY_ACID_METABOLISM", "REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION",
                       "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES","GO_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY", 
                       "GO_REGULATION_OF_CELL_CYCLE_PROCESS", "GO_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION")

pathways_df_gc <- pathway_prep(pathways, pathways_selected)

#enrichment  
tf_l_df_gc_so <- enricheR("NO3MGC - SO3MGC = 0", T, c_pw_gc, regulons_gc, pathways_df_gc)
#Jaccard index
tf_l_ji_gc_so <- JI_calculatoR(pathways_df_gc, regulons_gc)
tf_df_all_gc_so <- tf_selectoR(tf_l_df_gc_so, tf_l_ji_gc_so)
tf_ac <- read.csv("auc_gc", sep="\t")
tf_ac$Cell <- gsub("_rep", "", tf_ac$Cell)
rownames(tf_ac) <- tf_ac$Cell
tf_ac$Cell <- substr(tf_ac$Cell, 1, nchar(tf_ac$Cell)-4)
colnames(tf_ac)[2:ncol(tf_ac)]<- substr(colnames(tf_ac)[2:ncol(tf_ac)], 1, nchar(colnames(tf_ac)[2:ncol(tf_ac)])-3)

#SO vs NO
a <- data.frame("Cell"=names(clusters_tf), "an"=clusters_tf)
tf_ac_so <- merge(tf_ac, a, by.x=0,by.y="Cell")
cells <- tf_ac_so$an
tf_rn <- tf_ac_so$Row.names
tf_ac_so <- tf_ac_so[,colnames(tf_ac_so)%in%tf_df_all_gc_so$tf,]
rownames(tf_ac_so) <- tf_rn
tf = HeatmapAnnotation(group = cells,
                        col = list(group = c("YN" = "black", 
                                             "YS_good"=brewer.pal(9, "Greys")[6],
                                             "YS_bad"=brewer.pal(9, "Greys")[3])
                                   ))  
col <- colorRamp2(c(-3.5, 0, 4.5), c(brewer.pal(11, "RdBu")[10],"white", brewer.pal(11, "RdBu")[2]))
tf_ac_so_s <- scale(tf_ac_so)
tf_ac_so_s <- tf_ac_so_s[,order(factor(colnames(tf_ac_so_s),levels = tfs[38:1]))]
Heatmap(t(tf_ac_so_s),  col=col, show_column_names = T, top_annotation =tf, cluster_rows = F,show_row_names = T,
        heatmap_legend_param = list(title = ""))
tfs <- colnames(tf_ac_so_s)
save(tfs, file="tfs.Rdata")
