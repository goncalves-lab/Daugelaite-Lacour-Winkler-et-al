#Projection on superovulation PCA

library(Seurat)

load("ovulation.Rdata") #seurat object created by create_seuart_age_ov.R


# Oocytes

dge_oc_ov = read.csv("dge_oc_ov.csv") #output of dge.R

OC <- subset(ovulation, cell == "OC")
OC = SCTransform(OC)

#selecting genes that are significantly disregulated by ageing
counts = OC@assays$SCT@counts
counts_ov = counts[dge_oc_ov$Gene_id[dge_oc_ov$padj_y<0.05 & !is.na(dge_oc_ov$padj_y)],]

#computing PCA on young natural and superovulated OC
pca_ov = prcomp(t(as.matrix(counts_ov[, grep("3M", colnames(counts_ov))])), center = T, scale. = T)

#projecting old cells on superovulation PCA
pca_no = scale(t(as.matrix(counts_ov[, grep("NO12M", colnames(counts_ov))])), pca_ov$center, pca_ov$scale) %*% pca_ov$rotation
pca_so = scale(t(as.matrix(counts_ov[, grep("SO12M", colnames(counts_ov))])), pca_ov$center, pca_ov$scale) %*% pca_ov$rotation


#preparing data for plotting
df_plot_ov = data.frame(pc1 = pca_ov$x[,1], group = c("NO" = "NY", "SO" = "SY")[substr(rownames(pca_ov$x), 1, 2)])
df_plot_ov = rbind(df_plot_ov, data.frame(pc1 = pca_no[,1], group = "NO"))
df_plot_ov = rbind(df_plot_ov, data.frame(pc1 = pca_so[,1], group = "SO"))

#plots for figures
pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pc1_ovulation_stripchart.pdf", width = 5, height = 3)
stripchart(df_plot_ov$pc1 ~ factor(df_plot_ov$group, levels = unique(df_plot_ov$group)), 
           col = c("NY"= "black","NO"= "firebrick3", "SY"= "darkgrey", "SO"= "firebrick1")[unique(df_plot_ov$group)],
           pch = 20, las = 1, xlab = "superovulation PC1", method = "jitter")
dev.off()

means = aggregate(df_plot_ov$pc1, by = list(df_plot_ov$group), mean)

df_shifts = data.frame(group = c("SY - NY", "SO - NO", "SO - SY", "NO - NY"),
                       shift = -c(means$x[4] - means$x[2], means$x[3] - means$x[1],
                                  means$x[3] - means$x[4], means$x[1] - means$x[2]))[4:1,]
pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/shifts_pc1_ovulation.pdf", width = 5, height = 4)
stripchart(df_shifts$shift ~ factor(df_shifts$group, levels = df_shifts$group), pch = 19, las = 1, xlab = "superovulation PC1 distance", cex = 2)
abline(v = 0, lty = "dashed", col = "grey")
segments(x0 = c(df_shifts$shift), x1 = c(0,0,0,0), y0 = 1:4, y1 = 1:4)
arrows(x0 = df_shifts$shift[c(1,3)], x1 = df_shifts$shift[c(2,4)], y0 = c(1.5, 3.5), y1 = c(1.5, 3.5), col = "cyan4", code = 3, length = 0.1)
text(x = c(mean(df_shifts$shift[1:2]), mean(df_shifts$shift[3:4])), y = c(1.5, 3.5), pos = 3,
     labels = c("Interaction age x ovulation"), col = "cyan4")
dev.off()

pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pca_ovulation_background.pdf", width = 5, height = 4)
pc_percent = ((pca_ov$sdev^2)/sum(pca_ov$sdev^2))
plot(pca_ov$x[,1], pca_ov$x[,2], col = c("black", "darkgrey")[1+grepl("SO", rownames(pca_ov$x))], 
     pch = 19, xlab = paste0("PC1 (", round(100*pc_percent[1],1), "%)"), 
     ylab = paste0("PC2 (", round(100*pc_percent[2],1), "%)"), main = "Projection onto SY/NY PCA")
dev.off()


pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pca_ovulation_projection.pdf", width = 5, height = 4)
pc_percent = ((pca_ov$sdev^2)/sum(pca_ov$sdev^2))
plot(pca_ov$x[,1], pca_ov$x[,2], col = scales::alpha(c("black", "darkgrey")[1+grepl("SO", rownames(pca_ov$x))], 0.5), 
     pch = 19, xlab = paste0("PC1 (", round(100*pc_percent[1],1), "%)"), 
     ylab = paste0("PC2 (", round(100*pc_percent[2],1), "%)"), main = "Projection onto SY/NY PCA")
points(pca_no[,1], pca_no[,2], col = "firebrick3", pch = 19)
points(pca_so[,1], pca_so[,2], col = "firebrick1", pch = 19)
L = legend("topright", legend = c("", "", "Y", "O"), pch = 20, xpd = T,
           col = c("black","firebrick3", "darkgrey", "firebrick1"), bty = "n", ncol=2, inset = c(0,-0.15))
text(x = c(L$text$x[1]-0.15*L$rect$w, L$text$x[3]-0.15*L$rect$w), y=L$rect$top, label = c("N", "S"), xpd= NA)
dev.off()


##doing the same on GC

dge_gc_ov = read.csv("dge_gc_ov.csv") #output of dge.R

GC <- subset(ovulation, cell == "GC")
# GC = SCTransform(GC)

#selecting genes that are significantly disregulated by ageing
counts = GC@assays$SCT@counts
counts_ov = counts[rownames(counts)%in%dge_gc_ov$Gene_id[dge_gc_ov$padj_y<0.05 & !is.na(dge_gc_ov$padj_y)],]

#computing PCA on young natural and superovulated OC
pca_ov = prcomp(t(as.matrix(counts_ov[, grep("3M", colnames(counts_ov))])), center = T, scale. = T)

#projecting old cells on superovulation PCA
pca_no = scale(t(as.matrix(counts_ov[, grep("NO12M", colnames(counts_ov))])), pca_ov$center, pca_ov$scale) %*% pca_ov$rotation
pca_so = scale(t(as.matrix(counts_ov[, grep("SO12M", colnames(counts_ov))])), pca_ov$center, pca_ov$scale) %*% pca_ov$rotation

#preparing data for plotting
df_plot_ov = data.frame(pc1 = pca_ov$x[,1], group = c("NO" = "NY", "SO" = "SY")[substr(rownames(pca_ov$x), 1, 2)])
df_plot_ov = rbind(df_plot_ov, data.frame(pc1 = pca_no[,1], group = "NO"))
df_plot_ov = rbind(df_plot_ov, data.frame(pc1 = pca_so[,1], group = "SO"))

# plots
pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pc1_ovulation_stripchart_gc.pdf", width = 5, height = 3)
stripchart(df_plot_ov$pc1 ~ factor(df_plot_ov$group, levels = unique(df_plot_ov$group)), 
           col = c("NY"= "black","NO"= "firebrick3", "SY"= "darkgrey", "SO"= "firebrick1")[unique(df_plot_ov$group)],
           pch = 20, las = 1, xlab = "superovulation PC1", method = "jitter")
dev.off()

means = aggregate(df_plot_ov$pc1, by = list(df_plot_ov$group), mean)

df_shifts = data.frame(group = c("SY - NY", "SO - NO", "SO - SY", "NO - NY"),
                       shift = -c(means$x[4] - means$x[2], means$x[3] - means$x[1],
                                  means$x[3] - means$x[4], means$x[1] - means$x[2]))[4:1,]
pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/shifts_pc1_ovulation_gc.pdf", width = 5, height = 4)
stripchart(df_shifts$shift ~ factor(df_shifts$group, levels = df_shifts$group), pch = 19, las = 1, xlab = "superovulation PC1 distance", cex = 2)
abline(v = 0, lty = "dashed", col = "grey")
segments(x0 = c(df_shifts$shift), x1 = c(0,0,0,0), y0 = 1:4, y1 = 1:4)
arrows(x0 = df_shifts$shift[c(1,3)], x1 = df_shifts$shift[c(2,4)], y0 = c(1.5, 3.5), y1 = c(1.5, 3.5), col = "cyan4", code = 3, length = 0.1)
text(x = c(mean(df_shifts$shift[1:2]), mean(df_shifts$shift[3:4])), y = c(1.5, 3.5), pos = 3,
     labels = c("Interaction age x ovulation"), col = "cyan4")
dev.off()

pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pca_ovulation_background_gc.pdf", width = 5, height = 4)
pc_percent = ((pca_ov$sdev^2)/sum(pca_ov$sdev^2))
plot(pca_ov$x[,1], pca_ov$x[,2], col = c("black", "darkgrey")[1+grepl("SO", rownames(pca_ov$x))], 
     pch = 19, xlab = paste0("PC1 (", round(100*pc_percent[1],1), "%)"), 
     ylab = paste0("PC2 (", round(100*pc_percent[2],1), "%)"), main = "Projection onto SY/NY PCA")
dev.off()


pdf("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/figures/pca_ovulation_projection_gc.pdf", width = 5, height = 4)
pc_percent = ((pca_ov$sdev^2)/sum(pca_ov$sdev^2))
plot(pca_ov$x[,1], pca_ov$x[,2], col = scales::alpha(c("black", "darkgrey")[1+grepl("SO", rownames(pca_ov$x))], 0.5), 
     pch = 19, xlab = paste0("PC1 (", round(100*pc_percent[1],1), "%)"), 
     ylab = paste0("PC2 (", round(100*pc_percent[2],1), "%)"), main = "Projection onto SY/NY PCA")
points(pca_no[,1], pca_no[,2], col = "firebrick3", pch = 19)
points(pca_so[,1], pca_so[,2], col = "firebrick1", pch = 19)
L = legend("topright", legend = c("", "", "Y", "O"), pch = 20, xpd = T,
           col = c("black","firebrick3", "darkgrey", "firebrick1"), bty = "n", ncol=2, inset = c(0,-0.15))
text(x = c(L$text$x[1]-0.15*L$rect$w, L$text$x[3]-0.15*L$rect$w), y=L$rect$top, label = c("N", "S"), xpd= NA)
dev.off()


