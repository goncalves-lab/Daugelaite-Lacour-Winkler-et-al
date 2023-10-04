#this script calculates dShE in OC and GC
library(Seurat)
library(EntropyExplorer)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(data.table)
library(stringr)
library(biomaRt)
library(VennDiagram)
library(dplyr)
library(vioplot)


#functions
#####################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


data_prepareR <- function(counts, cell_type, group_name){
  counts_oc <- as.data.frame(counts[,grep(cell_type, colnames(counts))])
  if (cell_type=="OC"){
    names(clusters_consensus) <- gsub("GC", "OC", names(clusters_consensus))
  }
  #ny <- names(clusters_consensus[clusters_consensus==group_name])
  colnames(counts_oc) <- gsub("_rep", "", colnames(counts_oc))
  counts_oc_no_y <- counts_oc[,grep(group_name, colnames(counts_oc))]
  return(counts_oc_no_y)
}

EntropeR <- function(ov_df_1, ov_df_2){
  ov_ent <- EntropyExplorer(ov_df_1, ov_df_2, "dse", "bu", shift = c("auto","auto"))
  ov_ent <- as.data.frame(ov_ent)
  return(ov_ent)
}

plotteR <- function(ov_ent_fdr, group_1, group_2, sig){
  if (sig==T){
    ov_ent_fdr <- ov_ent_fdr[(ov_ent_fdr$`fdr p-value`<0.1),]} 
  a<- ggplot(ov_ent_fdr, aes(x=`SE(expm1)`, y=`SE(expm2)`)) + 
          geom_point(size=0.5)+
          theme_minimal()+
          theme(axis.text=element_text(size=20), axis.title=element_text(size=25))+
          labs(x=paste("ShE",group_1), y=paste("ShE",group_2))+
          theme(legend.position="right", 
                legend.title = element_text(colour="Black", size=18),
                legend.text = element_text(colour="black", size = 20))+
          geom_hline(yintercept=0, linetype="dashed", color = "red")+
          guides(color=guide_legend(override.aes=list(fill=NA),title=""))+
          geom_vline(xintercept=0, linetype="dashed", color = "red")+ geom_abline(intercept = 0, slope = 1)+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "white"))
  # Marginal density plot of x (top panel)
  group_1_df <- data.frame("ShE"=ov_ent_fdr$`SE(expm1)`, "group"=group_1)
  group_2_df <- data.frame("ShE"=ov_ent_fdr$`SE(expm2)`, "group"=group_2)
  group_df <- rbind(group_1_df, group_2_df)
  
  b <-ggplot(group_df, aes(ShE, fill=group)) + 
          geom_density(alpha=.5) + 
          scale_fill_manual(values = c('#999999','#E69F00')) + 
          theme(legend.position = "none")+
          theme_classic()+
          theme(axis.text=element_text(size=20), axis.title=element_text(size=25))+
          theme(legend.position="right", 
                legend.title = element_text(colour="Black", size=18),
                legend.text = element_text(colour="black", size = 20))
  print(a+b)
  
  # Marginal density plot of y (right panel)
}

#start
########
counts <- loadRData("NormalizedData.Rdata") #data were produced using Scnorm.R script

counts_oc_no_y <- data_prepareR(counts, "OC", "NO3M")
counts_oc_so <- data_prepareR(counts, "OC", "SO3M")
counts_oc_o <- data_prepareR(counts, "OC", "NO12M")

counts_gc_no_y <- data_prepareR(counts, "GC", "NO3M")
counts_gc_so <- data_prepareR(counts, "GC", "SO3M")
counts_gc_o <- data_prepareR(counts, "GC", "NO12M")

ent_y_no_so_oc <- EntropeR(counts_oc_so, counts_oc_no_y)
ent_y_no_o_oc <- EntropeR(counts_oc_o, counts_oc_no_y)


ent_y_no_so_gc <- EntropeR(counts_gc_so, counts_gc_no_y)
ent_y_no_o_gc <- EntropeR(counts_gc_o, counts_gc_no_y)

