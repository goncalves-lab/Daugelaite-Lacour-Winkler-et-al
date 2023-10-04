#performs permutation test to assess if regulons are differentially regulated
library(Seurat)
library(gsubfn)
library(dplyr)
library(coin)
library("FSA")
library(rcompanion)
library(multcompView)
library(rcompanion)
library(superheat)
library(RColorBrewer)
library(purrr)
library(tidyverse)



#permutation of cell types
tissue_permutoR <- function(metadata, scenic_df, cell_type){
  metadata_t <- metadata[metadata$cell_type==cell_type,]
  scenic_df_t <- scenic_df[rownames(scenic_df)%in%metadata_t$cell_name,]
  scenic_df_t <- merge(metadata_t, scenic_df_t, by.x="cell_name", by.y=0)
  test_vec <- c()
  #PT_total <- data.frame("Comparison"="counter", "Stat"=1, "p.value"=1, "p.adjust"=1,"tf"="counter") 
  for (i in seq(from=4, to=ncol(scenic_df_t))){
    tf <- as.data.frame(scenic_df_t[,c(1,2,i)])
    rownames(tf) <- scenic_df_t$cell_name
    colnames(tf)[3] <- "tf_name"
    tf$condition <- as.factor(tf$condition)
    if (table(tf$tf_name)[[1]]<(0.9*nrow(tf))){
      test1 <- independence_test( tf_name ~ condition, data=tf, distribution=approximate(10000))
      test_vec_t <- pvalue(test1)[1]
      names(test_vec_t) <- colnames(scenic_df_t)[i]
      test_vec <- c(test_vec, test_vec_t)}
    #PT <- pairwisePermutationTest(tf_name ~ cell_types, data=tf,
    #                             method="fdr")
    #PT$tf <- colnames(scenic_df_m)[i]
    #PT_total <- rbind(PT_total,PT)
  }
  #test_vec <- test_vec[test_vec<0.05]
  return(test_vec)
}

tissue_permutoR_pairwise <- function(metadata, scenic_df, cell_type, sig_scenic){
  metadata_t <- metadata[metadata$cell_type==cell_type,]
  scenic_df_t <- scenic_df[rownames(scenic_df)%in%metadata_t$cell_name,]
  scenic_df_t <- scenic_df_t[,colnames(scenic_df_t)%in%sig_scenic$tf]
  scenic_df_t <- merge(metadata_t, scenic_df_t, by.x="cell_name", by.y=0)
  
  test_vec <- c()
  PT_total <- data.frame("Comparison"="counter", "Stat"=1, "p.value"=1, "p.adjust"=1,"tf"="counter") 
  for (i in seq(from=4, to=ncol(scenic_df_t))){
    tf <- as.data.frame(scenic_df_t[,c(1,2,i)])
    rownames(tf) <- rownames(scenic_df_t)
    colnames(tf)[3] <- "tf_name"
    tf$condition <- as.factor(tf$condition)
    if (table(tf$tf_name)[[1]]<(0.9*nrow(tf))){
      #test1 <- independence_test( tf_name ~ phase + sample, data=tf, distribution=approximate(10000))
      #test_vec_t <- pvalue(test1)[1]
      #names(test_vec_t) <- colnames(scenic_df_t)[i]
      #test_vec <- c(test_vec, test_vec_t)
      PT <- pairwisePermutationTest(tf_name ~ condition, data=tf,
                                    method="fdr")
      PT$tf <- colnames(scenic_df_t)[i]
      PT_total <- rbind(PT_total,PT)}
  }
  #test_vec <- test_vec[test_vec<0.05]
  return(PT_total)
}


##############
scenic_df <- read.csv("auc_gc", header=T, sep="", row.names = 1) #produced in script scenic.R
rownames(scenic_df) <- gsub("S0", "SO", rownames(scenic_df))
metadata <- data.frame(cell_name=rownames(scenic_df), condition=gsub("_rep", "",rownames(scenic_df)))
metadata$condition <- substr(metadata$condition, 1, nchar(metadata$condition)-2)
metadata$cell_type <- substr(metadata$condition, nchar(metadata$condition)-1, nchar(metadata$condition))

#OC  
oc_result <- tissue_permutoR(metadata, scenic_df, "OC")
oc_result_df <- data.frame("tf"= names(oc_result), "p_value" = unname(oc_result), "p_value_adj" = p.adjust(unname(oc_result), method="BH"))
oc_result_df_sig <- oc_result_df[oc_result_df$p_value_adj<0.05,]
oc_pw <- tissue_permutoR_pairwise(metadata, scenic_df, "OC", oc_result_df_sig)
oc_pw <- oc_pw[2:nrow(oc_pw),]
oc_pw <- oc_pw[, c(1,4,5)]
oc_pw <- oc_pw[,c(1,3,2)]
colnames(oc_pw)[3] <- "p_value"
oc_pw$p_adj <- p.adjust(oc_pw$p_value, method="BH")

#GC  
gc_result <- tissue_permutoR(metadata, scenic_df, "GC")
gc_result_df <- data.frame("tf"= names(gc_result), "p_value" = unname(gc_result), "p_value_adj" = p.adjust(unname(gc_result), method="BH"))
gc_result_df_sig <- gc_result_df[gc_result_df$p_value_adj<0.05,]
gc_pw <- tissue_permutoR_pairwise(metadata, scenic_df, "GC", gc_result_df_sig)
gc_pw <- gc_pw[2:nrow(gc_pw),]
gc_pw <- gc_pw[, c(1,4,5)]
gc_pw <- gc_pw[,c(1,3,2)]
colnames(gc_pw)[3] <- "p_value"
gc_pw$p_adj <- p.adjust(gc_pw$p_value, method="BH")

save(oc_pw, file="oc_separate_pw_nl.Rdata")
save(gc_pw, file="gc_separate_pw_nl.Rdata")
