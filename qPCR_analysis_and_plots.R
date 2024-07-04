# Analysis of qPCR values

#-----Libraries-------

library(dplyr)
library(tidyr)

#------Data--------

qpcr = read.csv("qPCR_values.csv", sep = ";", dec = ",") # raw values from qPCR machine, each row is a well
metadata = read.csv("qPCR_samples.csv", sep = ";", dec = ",") #metadata table

#removing technical replicates (the average between both is computed by the qPCR machine)
qpcr = qpcr[!duplicated(qpcr[,c(1,4,5)]),]

#normalising value relative to GAPDH
dCt = qpcr %>% 
  group_by(Plate, Sample.Name) %>%
  mutate(dCt = Ct.Mean - Ct.Mean[Target.Name == "GAPDH"]) %>%
  filter(Ct.SD < 0.25)

#adding metadata information
dCt$Sample = metadata$Cell[match(paste0(dCt$Plate, dCt$Sample.Name), paste0(metadata$Plate, metadata$Sample))]
dCt = dCt[, c("Sample", "Target.Name", "dCt")]

#changing format of the table to ease representation
gene_ct = dCt %>% spread(key = Sample, value = dCt) %>% as.data.frame()
rownames(gene_ct) = gene_ct$Target.Name
gene_ct = as.data.frame(t(gene_ct[,-1]))

gene_values = 2**-gene_ct
gene_values$condition = substr(rownames(gene_values), 1,2)
gene_values = gene_values[, colnames(gene_values)!="GAPDH"]


#adding 
load("gc_class_all.Rdata") #classification information from the classifier scripts (the same samples were sequenced using Smart-Seq2)
gene_values$group = gc_class_all[tolower(rownames(gene_values))]
gene_values$group[gene_values$condition=="NO"] = "NY"
gene_ct$group = gc_class_all[tolower(rownames(gene_ct))]
gene_ct$group[gene_ct$condition=="NO"] = "NY"
gene_ct$condition = substr(rownames(gene_ct), 1,2)



#plots for the figures
pdf("qPCR_all.pdf", width = 10, height = 7)
bp = gene_values %>% filter(condition == "SO", group == "YS_bad") %>% 
  dplyr::select(!c(condition, group)) %>% 
  boxplot(pch = 20, log = "y", col = "#6ED2AC", at = 1:12-0.15, boxwex = 0.12, 
          xaxt = 'n', las = 1, ylim = c(min(gene_values[,1:12], na.rm = T), max(gene_values[,1:12], na.rm = T)*5))
gene_values %>% filter(condition == "SO", group == "YS_good") %>% 
  dplyr::select(!c(condition, group)) %>% 
  boxplot(pch = 20, log = "y", col = "#236B4F", at = 1:12, boxwex = 0.12, add = T, axes = F)
gene_values %>% filter(group == "NY") %>% 
  dplyr::select(!c(condition, group)) %>% 
  boxplot(pch = 20, log = "y", col = "grey20", at = 1:12+0.15, boxwex = 0.12, add = T, axes = F)
axis(1, at = 1:12, labels = F)
text(1:12, labels = colnames(gene_values)[1:12], srt = 45, y = 5e-8, xpd = NA, adj = 1)

pvals = apply(gene_values[,1:12], 2, function(x){wilcox.test(x[gene_values$condition=="SO" & gene_values$group!="NA"] ~ 
                                                               gene_values$group[gene_values$condition=="SO" & gene_values$group!="NA"])$p.value})
text(x = 1:12, y = 4*max(c(bp$stats, bp$out)), pos = 1, cex = 1,
     labels = c("***", "**", "*", "ns")[cut(pvals, breaks = c(0,0.001,0.01,0.05,1), labels = 1:4)], xpd = T)
dev.off()


pdf("qPCR_all_dots.pdf", width = 10, height = 7) #Supplementary Figure 3f
gene_values %>% filter(condition == "SO", group == "YS_bad") %>% 
  dplyr::select(!c(condition, group)) %>% 
  stripchart(pch = 20, log = "y", col = ggplot2::alpha("#6ED2AC", 0.5), at = 1:12+0.2, boxwex = 0.12, 
             vertical = T, method = "jitter", jitter = 0.02,
             xaxt = 'n', las = 1, ylim = c(min(gene_values[,1:12], na.rm = T), max(gene_values[,1:12], na.rm = T)*5))
gene_values %>% filter(condition == "SO", group == "YS_good") %>% 
  dplyr::select(!c(condition, group)) %>% 
  stripchart(pch = 20, log = "y", col = ggplot2::alpha("#236B4F", 0.5), at = 1:12, boxwex = 0.12, 
             add = T, axes = F, 
             vertical = T, method = "jitter", jitter = 0.02)
gene_values %>% filter(group == "NY") %>% 
  dplyr::select(!c(condition, group)) %>% 
  stripchart(pch = 20, log = "y", col = ggplot2::alpha("black", 0.5), at = 1:12-0.2, boxwex = 0.12, 
             add = T, axes = F, 
             vertical = T, method = "jitter", jitter = 0.02)
axis(1, at = 1:12, labels = F)
text(1:12, labels = colnames(gene_values)[1:12], srt = 45, y = 5e-8, xpd = NA, adj = 1)

meds = gene_values %>% dplyr::select(!c(condition)) %>% 
  group_by(group) %>% summarise(across(.fns = median, na.rm=T))
meds = meds[match(c("YS_bad", "YS_good", "NY"), meds$group),]
points(x = 1:12+0.2, pch = 23,
       y = as.numeric(meds[meds$group=="YS_bad", 2:ncol(meds)]), 
       col = "#6ED2AC", lwd = 3, bg = "white")
points(x = 1:12, pch = 23,
       y = as.numeric(meds[meds$group=="YS_good", 2:ncol(meds)]), 
       col = "#236B4F", lwd = 3, bg = "white")
points(x = 1:12-0.2, pch = 23,
       y = as.numeric(meds[meds$group=="NY", 2:ncol(meds)]), 
       col = "black", lwd = 3, bg = "white")

pvals = apply(gene_values[,1:12], 2, function(x){wilcox.test(x[gene_values$condition=="SO" & gene_values$group!="NA"] ~ 
                                                               gene_values$group[gene_values$condition=="SO" & gene_values$group!="NA"])$p.value})
text(x = 1:12, y = 4*max(gene_values[,1:12], na.rm = T), pos = 1, cex = 1,
     labels = c("***", "**", "*", "ns")[cut(pvals, breaks = c(0,0.001,0.01,0.05,1), labels = 1:4)], xpd = T)
dev.off()


pdf("qPCR_subset.pdf", width = 10, height = 7) #Main Figure 4e
genes_oi = names(pvals)[pvals<0.05]
n = length(genes_oi)
bp = gene_values %>% filter(condition == "SO", group == "YS_bad") %>% 
  dplyr::select(!c(condition, group)) %>% 
  dplyr::select(all_of(genes_oi)) %>%
  boxplot(pch = 20, log = "y", col = "#6ED2AC", at = 1:n+0.15, boxwex = 0.12, 
          xaxt = 'n', las = 1, ylim = c(min(gene_values[,1:12], na.rm = T), max(gene_values[,1:12], na.rm = T)*5))
gene_values %>% filter(condition == "SO", group == "YS_good") %>% 
  dplyr::select(!c(condition, group)) %>% 
  dplyr::select(all_of(genes_oi)) %>%
  boxplot(pch = 20, log = "y", col = "#236B4F", at = 1:n, boxwex = 0.12, add = T, axes = F)
gene_values %>% filter(group == "NY") %>% 
  dplyr::select(!c(condition, group)) %>% 
  dplyr::select(all_of(genes_oi)) %>%
  boxplot(pch = 20, log = "y", col = "grey20", at = 1:n-0.15, boxwex = 0.12, add = T, axes = F)
axis(1, at = 1:n, labels = F)
text(1:n, labels = genes_oi, srt = 45, y = 5e-8, xpd = NA, adj = 1)
text(x = 1:n, y = 4*max(c(bp$stats, bp$out)), pos = 1, cex = 2,
     labels = c("***", "**", "*", "ns")[cut(pvals[genes_oi], breaks = c(0,0.001,0.01,0.05,1), labels = 1:4)], xpd = T)
dev.off()

