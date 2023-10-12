#Script for cell-to-cell communication between oocytes and granulosa cells

#-----Library--------
library(Seurat)
library(biomaRt)
library(hash)


#-----Parameters and data--------

ensembl_mouse=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list_mouse=getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=ensembl_mouse)



#Preparing the ligand-receptor data
LR_pairs_CellChat = read.csv("Mouse-2020-Jin-LR-pairs.csv") #CellChat DB downloaded from their website on June 22nd 2021, also found as CellChatDB.mouse in the CellChat package

## We parse the interaction name to get the gene names for the ligand and the receptor genes
pairs_CellChat = LR_pairs_CellChat$interaction_name_2
pairs_CellChat = gsub(" +"," ",pairs_CellChat)
pairs_CellChat = strsplit(pairs_CellChat," - ")
pairs_CellChat = lapply(pairs_CellChat, function(X){X[2]=gsub("(","",X[2], fixed = T); return(X)})
pairs_CellChat = lapply(pairs_CellChat, function(X){X[2]=gsub(")","",X[2], fixed = T); return(X)})
pairs_CellChat = lapply(pairs_CellChat, function(X){X[2]=strsplit(X[2],"+", fixed = T); return(X)})

LR_df_CellChat = data.frame(ligand = sapply(pairs_CellChat, function(X){return(X[[1]])}),
                            receptor1 = sapply(pairs_CellChat, function(X){return(X[[2]][1])}),
                            receptor2 = sapply(pairs_CellChat, function(X){return(X[[2]][2])}))


#Preparing the pairing metadata
load("ovulation.Rdata") #Seurat object created by script create_seurat_age_ov.R

paired_id = data.frame() # This table will contain the list of paired oocytes and granulosa cells alongside their age / treatment group
gc_list = names(ovulation$label)[grep("GC",names(ovulation$label))]
oc_list = names(ovulation$label)[grep("OC",names(ovulation$label))]
for (i in 1:length(gc_list)) {
  if (length(grep(gsub("GC","OC",gc_list[i]),oc_list))) {
    paired_id = rbind(paired_id, c(gc_list[i], oc_list[grep(gsub("GC","OC",gc_list[i]),oc_list)]))
  }
}
colnames(paired_id) = c("GC", "OC")
paired_id$group = paste(ovulation$age[paired_id$GC],ovulation$ovulation[paired_id$GC],sep="_")



#Preparing the count data
load("NormalizedData.Rdata") #count data normalized accounting for gene length by SCnorm
count_matrix = NormalizedData
rownames(count_matrix) = ensemnl_list_mouse$external_gene_name[match(rownames(count_matrix), ensemnl_list_mouse$ensembl_gene_id)]
count_matrix = count_matrix[,c(paired_id$OC,paired_id$GC)]
count_matrix_GC = count_matrix[,paired_id$GC]
count_matrix_OC = count_matrix[,paired_id$OC]

## Filtering for genes that are expressed in at least 10 cells of each cell type
count_matrix_GC = count_matrix_GC[rowSums(count_matrix_GC>0)>10,]
count_matrix_OC = count_matrix_OC[rowSums(count_matrix_OC>0)>10,]


ligand_matrix_GC = count_matrix_GC[rownames(count_matrix_GC)%in%LR_df_CellChat$ligand,]
receptor_matrix_GC = count_matrix_GC[rownames(count_matrix_GC)%in%c(LR_df_CellChat$receptor1, LR_df_CellChat$receptor2),]
ligand_matrix_OC = count_matrix_OC[rownames(count_matrix_OC)%in%LR_df_CellChat$ligand,]
receptor_matrix_OC = count_matrix_OC[rownames(count_matrix_OC)%in%c(LR_df_CellChat$receptor1, LR_df_CellChat$receptor2),]

# Filtering the tables to keep the top 100 most expressed ligand and receptor genes in each cell type
ligand_matrix_GC = ligand_matrix_GC[order(-rowMeans(ligand_matrix_GC, na.rm = F))[1:100],]
receptor_matrix_GC = receptor_matrix_GC[order(-rowMeans(receptor_matrix_GC, na.rm = F))[1:100],]
ligand_matrix_OC = ligand_matrix_OC[order(-rowMeans(ligand_matrix_OC, na.rm = F))[1:100],]
receptor_matrix_OC = receptor_matrix_OC[order(-rowMeans(receptor_matrix_OC, na.rm = F))[1:100],]

count_matrix_GC = rbind(ligand_matrix_GC, receptor_matrix_GC)
count_matrix_OC = rbind(ligand_matrix_OC, receptor_matrix_OC)

#-----Main-------

df_LR_CellChat = matrix(ncol=2*length(paired_id$GC), nrow=nrow(LR_pairs_CellChat))
df_autocrine_CellChat = data.frame()
for (i in 1:nrow(LR_pairs_CellChat)) {
  # Parsing the interaction name to get the ligand and receptor gene names
  pair = LR_pairs_CellChat$interaction_name_2[i]
  pair = gsub(" +"," ",pair)
  pair = strsplit(pair," - ")
  pair[[1]][2] = gsub("(","",pair[[1]][2], fixed = T)
  pair[[1]][2] = gsub(")","",pair[[1]][2], fixed = T)
  pair[[1]][2] = strsplit(pair[[1]][2],"+", fixed = T)
  ligand = pair[[1]][[1]]
  receptor = pair[[1]][[2]]
  if (receptor=="Gp complex") {receptor = c("Gp1ba", "Gp1bb", "Gp5", "Gp9")}
  
  # Extracting the counts from the count table for each cell type
  ligand_score_OC = count_matrix_OC[rownames(count_matrix_OC)%in%ligand,]
  receptor_score_OC = count_matrix_OC[rownames(count_matrix_OC)%in%receptor,]
  ligand_score_GC = count_matrix_GC[rownames(count_matrix_GC)%in%ligand,]
  receptor_score_GC = count_matrix_GC[rownames(count_matrix_GC)%in%receptor,]
  
  if (length(receptor) >1) { # if there are two subunits of the receptor
    #if at least one subunit is not expressed, receptor score is NA
    if (is.null(nrow(receptor_score_OC)) || nrow(receptor_score_OC)!=length(receptor)) 
    {receptor_score_OC = rep(NA,nrow(paired_id))} 
    # if both subunits are expressed, the least expressed is considered (it is the limiting factor assuming 1:1 ratio)
    else {receptor_score_OC = apply(receptor_score_OC, MARGIN = 2, min)}
    # we do the same in GC
    if (is.null(nrow(receptor_score_GC)) || nrow(receptor_score_GC)!=length(receptor))
      {receptor_score_GC = rep(NA,nrow(paired_id))}
    else {receptor_score_GC = apply(receptor_score_GC, MARGIN = 2, min)}
    
  } else { #if there is only one receptor
    # if the receptor is not expressed, its score is NA
    if (!is.null(nrow(receptor_score_OC))) {receptor_score_OC = rep(NA,nrow(paired_id))}
    if (!is.null(nrow(receptor_score_GC))) {receptor_score_GC = rep(NA,nrow(paired_id))}
  }
  
  # if the ligand is not expressed, its score is NA
  if (!is.null(nrow(ligand_score_GC))) {ligand_score_GC = rep(NA,nrow(paired_id))}
  if (!is.null(nrow(ligand_score_OC))) {ligand_score_OC = rep(NA,nrow(paired_id))}
  
  # the interaction score is the product of the ligand and receptor scores
  newline = c(ligand_score_OC*receptor_score_GC, receptor_score_OC*ligand_score_GC)
  df_LR_CellChat[i,] = newline
  df_autocrine_CellChat = rbind(df_autocrine_CellChat, c(ligand_score_OC*receptor_score_OC, receptor_score_GC*ligand_score_GC))
}

df_LR_CellChat = as.data.frame(df_LR_CellChat)
colnames(df_LR_CellChat) = c(paired_id$OC,paired_id$GC) #The colnames are the labels of the cells producing the ligand
rownames(df_LR_CellChat) = LR_pairs_CellChat$interaction_name

colnames(df_autocrine_CellChat) = c(paired_id$OC,paired_id$GC)
rownames(df_autocrine_CellChat) = LR_pairs_CellChat$interaction_name



## testing difference between groups of age / treatment using a wilcoxon test
df_pvalues_CellChat = data.frame()
for (i in 1:nrow(df_LR_CellChat)) {
  # if the interaction score is NA, the pvalue is also NA
  if (sum(!is.na(df_LR_CellChat[i,paired_id$GC]))==0) {GC_OC_SO_Y_O = GC_OC_NO_Y_O = GC_OC_Y_SO_NO = GC_OC_O_SO_NO = NA}
  else {
    GC_OC_SO_Y_O = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="O_SO"]])), 
                               y = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="Y_SO"]])))$p.value
    GC_OC_NO_Y_O = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="O_NO"]])), 
                               y = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="Y_NO"]])))$p.value
    GC_OC_Y_SO_NO = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="Y_SO"]])), 
                                y = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="Y_NO"]])))$p.value
    GC_OC_O_SO_NO = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="O_SO"]])), 
                                y = log(as.numeric(df_LR_CellChat[i,paired_id$GC[paired_id$group=="O_NO"]])))$p.value
  }
  
  if (sum(!is.na(df_LR_CellChat[i,paired_id$OC]))==0) {OC_GC_SO_Y_O = OC_GC_NO_Y_O = OC_GC_Y_SO_NO = OC_GC_O_SO_NO = NA}
  else {
    OC_GC_SO_Y_O = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="O_SO"]])), 
                               y = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="Y_SO"]])))$p.value
    OC_GC_NO_Y_O = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="O_NO"]])), 
                               y = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="Y_NO"]])))$p.value
    OC_GC_Y_SO_NO = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="Y_SO"]])), 
                                y = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="Y_NO"]])))$p.value
    OC_GC_O_SO_NO = wilcox.test(x = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="O_SO"]])), 
                                y = log(as.numeric(df_LR_CellChat[i,paired_id$OC[paired_id$group=="O_NO"]])))$p.value
    
  }
  
  if (sum(!is.na(df_autocrine_CellChat[i,paired_id$OC]))==0) {OC_OC_SO_Y_O = OC_OC_NO_Y_O = OC_OC_Y_SO_NO = OC_OC_O_SO_NO = NA}
  else {
    OC_OC_SO_Y_O = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="O_SO"]])), 
                               y = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="Y_SO"]])))$p.value
    OC_OC_NO_Y_O = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="O_NO"]])), 
                               y = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="Y_NO"]])))$p.value
    OC_OC_Y_SO_NO = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="Y_SO"]])), 
                                y = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="Y_NO"]])))$p.value
    OC_OC_O_SO_NO = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="O_SO"]])), 
                                y = log(as.numeric(df_autocrine_CellChat[i,paired_id$OC[paired_id$group=="O_NO"]])))$p.value
    
  }
  if (sum(!is.na(df_autocrine_CellChat[i,paired_id$GC]))==0) {GC_GC_SO_Y_O = GC_GC_NO_Y_O = GC_GC_Y_SO_NO = GC_GC_O_SO_NO = NA}
  else {
    GC_GC_SO_Y_O = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="O_SO"]])), 
                               y = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="Y_SO"]])))$p.value
    GC_GC_NO_Y_O = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="O_NO"]])), 
                               y = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="Y_NO"]])))$p.value
    GC_GC_Y_SO_NO = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="Y_SO"]])), 
                                y = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="Y_NO"]])))$p.value
    GC_GC_O_SO_NO = wilcox.test(x = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="O_SO"]])), 
                                y = log(as.numeric(df_autocrine_CellChat[i,paired_id$GC[paired_id$group=="O_NO"]])))$p.value
    
  }
  
  df_pvalues_CellChat = rbind(df_pvalues_CellChat,c(OC_GC_Y_SO_NO,OC_GC_O_SO_NO,OC_GC_NO_Y_O,OC_GC_SO_Y_O,
                                                    GC_OC_Y_SO_NO,GC_OC_O_SO_NO,GC_OC_NO_Y_O,GC_OC_SO_Y_O,
                                                    OC_OC_Y_SO_NO,OC_OC_O_SO_NO,OC_OC_NO_Y_O,OC_OC_SO_Y_O,
                                                    GC_GC_Y_SO_NO,GC_GC_O_SO_NO,GC_GC_NO_Y_O,GC_GC_SO_Y_O))
}

rownames(df_pvalues_CellChat) = rownames(df_LR_CellChat)
colnames(df_pvalues_CellChat) = c("OC_GC_Y_SO_NO","OC_GC_O_SO_NO","OC_GC_NO_Y_O","OC_GC_SO_Y_O",
                                  "GC_OC_Y_SO_NO","GC_OC_O_SO_NO","GC_OC_NO_Y_O","GC_OC_SO_Y_O",
                                  "OC_OC_Y_SO_NO","OC_OC_O_SO_NO","OC_OC_NO_Y_O","OC_OC_SO_Y_O",
                                  "GC_GC_Y_SO_NO","GC_GC_O_SO_NO","GC_GC_NO_Y_O","GC_GC_SO_Y_O")

df_padj_CellChat = as.data.frame(matrix(p.adjust(as.matrix(df_pvalues_CellChat), method = "fdr"), ncol = dim(df_pvalues_CellChat)[2]))
rownames(df_padj_CellChat) = rownames(df_pvalues_CellChat)
colnames(df_padj_CellChat) = colnames(df_pvalues_CellChat)




# Saving results
write.table(df_LR_CellChat, "LR_expression_product_paracrine.tsv", col.names = T, row.names = T, quote = F)
write.table(df_autocrine_CellChat, "LR_expression_product_autocrine.tsv", col.names = T, row.names = T, quote = F)
write.table(df_padj_CellChat, "LR_padj_allconditions.tsv", col.names = T, row.names = T, quote = F)
