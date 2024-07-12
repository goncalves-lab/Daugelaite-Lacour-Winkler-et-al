#data prep for tf-based classifier
library(Seurat)
library(ggplot2)
library(caret)
library(glmmTMB)

############################
#TFs data preparation
#produced in scenic.R script
auc_df <- read.csv("auc_matrix_joined.csv", 
                   row.names = 1)
rownames(auc_df) <- gsub("\\(.*\\)", "", rownames(auc_df))
rownames(auc_df) <- gsub("_extended", "", rownames(auc_df))

load("tfs.Rdata")#produced in scenic_post.R script
auc_df <- as.data.frame(t(auc_df))
colnames(auc_df) <- gsub(" ", "", colnames(auc_df))
ori_rownames = ori_rownames[colnames(auc_df)%in%tfs]
auc_df <- auc_df[,colnames(auc_df)%in%tfs]

#identifying TF which have significantly different activity in SO good and bad
load("clusters_Y_consensus.Rdata") #consensus clusters
clusters_consensus <- clusters_consensus[!is.na(clusters_consensus)]
clusters_consensus <- clusters_consensus[grep("SO", names(clusters_consensus))]
clusters_consensus=droplevels(clusters_consensus)

auc_tf <- auc_df[rownames(auc_df)%in%names(clusters_consensus),]
clusters_df <- data.frame(cells=names(clusters_consensus), label=clusters_consensus)
auc_tf <- merge(auc_tf, clusters_df, by.x=0, by.y="cells")

load("ovulation.Rdata") #produced using scenic.R script
mouse_df <- data.frame(cells=colnames(ovulation), mouse=ovulation$mouse)
auc_tf <- merge(auc_tf, mouse_df, by.x="Row.names", by.y="cells")

auc_tf = auc_tf[, colSums(auc_tf==0)==0]

tf_df_final <- data.frame()
for(i in seq(from=2, to=ncol(auc_tf)-2)){
  tf_sig <- auc_tf[,c(i, ncol(auc_tf)-1, ncol(auc_tf))]
  tf_name <- colnames(auc_tf)[i]
  colnames(tf_sig)[1] <- "tf"
  tf_sig <- glmmTMB(tf~label+(1 | mouse), data = tf_sig, family=beta_family())
  tf_df <- data.frame("tf"=tf_name, "p_value"=summary(tf_sig)[6][[1]][[1]][8])
  tf_df_final <- rbind(tf_df, tf_df_final)
}

tf_df_final$p_adj <- p.adjust(tf_df_final$p_value, method = "fdr")
tf_df_final <- tf_df_final[tf_df_final$p_adj<0.05,]

ori_rownames = ori_rownames[colnames(auc_df)%in%tf_df_final$tf]
auc_df <- auc_df[,colnames(auc_df)%in%tf_df_final$tf]

#identifying near-zero variance predictors
nzv <- nearZeroVar(auc_df, saveMetrics= TRUE)

#identifying correlated predictors
descrCor <- cor(auc_df)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .8)
ori_rownames = ori_rownames[-highlyCorDescr]
auc_df <- auc_df[,-highlyCorDescr]

#identifying linear dependancies
comboInfo <- findLinearCombos(auc_df)

#splitting dataset

auc_df_so <- auc_df[rownames(auc_df)%in%clusters_df$cells,]
auc_df_so <- auc_df_so[match(clusters_df$cells, rownames(auc_df_so)),]
set.seed(3456)
trainIndex <- createDataPartition(clusters_df$label, p = .7, 
                                  list = FALSE, 
                                  times = 1)
auc_Train <- auc_df_so[ trainIndex,]
auc_Test  <- auc_df_so[-trainIndex,]

#center and scaling
preProcValues <- preProcess(auc_Train, method = c("center", "scale"))

auc_trainTransformed <- predict(preProcValues, auc_Train)
auc_testTransformed <- predict(preProcValues, auc_Test)


auc_trainTransformed <- merge(auc_trainTransformed, clusters_df, by.x=0, by.y="cells")
rownames(auc_trainTransformed) <- auc_trainTransformed$Row.names
auc_trainTransformed <- auc_trainTransformed[,2:ncol(auc_trainTransformed)]
colnames(auc_trainTransformed) = c(ori_rownames, "label")
save(auc_trainTransformed, file="auc_trainTransformed.Rdata")

auc_testTransformed <- merge(auc_testTransformed, clusters_df, by.x=0, by.y="cells")
rownames(auc_testTransformed) <- auc_testTransformed$Row.names
auc_testTransformed <- auc_testTransformed[,2:ncol(auc_testTransformed)]
colnames(auc_testTransformed) = c(ori_rownames, "label")

save(auc_testTransformed, file="auc_testTransformed.Rdata")




