# Analysis of HCR measurements
# This is done using the values measured using Fiji macro

#-----Libraries-------

library(dplyr)
library(lme4)


#--------Main--------

#for each COC, cells are segmented based on DAPI and the other channels are quantified
measure = read.csv("hcr_coc_measurements.csv", sep = ",", dec = ".") #values from Fiji
chan_vect = c("FAR-RED" = "Areg", "RED" = "Esr2", "GREEN" = "Npr2")

#parsing names for metadata information
measure$mouse = sapply(measure$Label, function(x){strsplit(x, "_")[[1]][1]})
measure$ovulation = sapply(measure$Label, function(x){strsplit(x, "_")[[1]][2]})
measure$coc = sapply(measure$Label, function(x){strsplit(x, "[_-]")[[1]][3]})
measure$coc = sapply(measure$coc, function(x){strsplit(x, "\\.")[[1]][1]})
measure$cell = measure$X
measure$channel = chan_vect[sapply(measure$Label, function(x){strsplit(x, " ")[[1]][2]})]

#removing cells that are more than 2 standard deviation from the mean
measure$cellID = paste(measure$mouse, measure$coc, measure$cell)
measure = measure %>%
  group_by(mouse, channel, coc) %>%
  mutate(outlier = Mean - mean(Mean, na.rm = T) > 2*sd(Mean, na.rm = T))

outlierID = unique(measure$cellID[measure$outlier])

good_cells = measure %>%
  filter(!cellID %in% outlierID)

#Averaging the cells in each COC (to be closer to the small bulks in Smart-seq2 experiments)
norm_per_coc = good_cells %>%
  group_by(mouse, channel, coc, ovulation) %>%
  summarise(norm_mean = mean(Mean, na.rm = T))
#Normalising the values based on the non-stained control COC
norm_per_coc = norm_per_coc %>%
  group_by(mouse, channel) %>%
  mutate(norm_mean = norm_mean - mean(norm_mean[coc=="no"], na.rm = T))




###fitting mixed models
mixed.Areg = norm_per_coc %>% filter(channel == "Areg", coc != "no") %>%
  lmer(norm_mean ~ ovulation + (1|mouse), data = ., REML = F) #singular value error

mixed.Npr2 = norm_per_coc %>% filter(channel == "Npr2", coc != "no") %>%
  lmer(norm_mean ~ ovulation + (1|mouse), data = ., REML = F)  #singular value error

mixed.Esr2 = norm_per_coc %>% filter(channel == "Esr2", coc != "no") %>%
  lmer(norm_mean ~ ovulation + (1|mouse), data = ., REML = F)
reduced.Esr2 = norm_per_coc %>% filter(channel == "Esr2", coc != "no") %>%
  lmer(norm_mean ~ 1 + (1|mouse), data = ., REML = F)
summary(mixed.Esr2)
anova(reduced.Esr2, mixed.Esr2)






### Figures for paper

mixed.Esr2 = norm_per_coc %>% filter(channel == "Esr2", coc != "no") %>%
  lmer(norm_mean ~ ovulation + (1|mouse), data = ., REML = F)
reduced.Esr2 = norm_per_coc %>% filter(channel == "Esr2", coc != "no") %>%
  lmer(norm_mean ~ 1 + (1|mouse), data = ., REML = F)
summary(mixed.Esr2)
aov = anova(reduced.Esr2, mixed.Esr2)

pdf("hcr_esr2_boxplot.pdf", width = 4, height = 5) #Supplementary Figure 1f
norm_per_coc %>% filter(channel=="Esr2", coc!="no") %>% 
  boxplot(norm_mean ~ mouse + channel, ., pch = 20, 
          las = 2, main = "Esr2", frame = F,
          col = ggplot2::alpha(c("black", "#727272", "black"),0.7), 
          at = c(0.8,2, 1.2), boxwex = 0.3, 
          lty = 1, staplelty = 0, whisklwd = 2, medlwd = 3,
          xlab = NA, xaxt = "n", ylab = c("Normalized mean intensity per COC"))
text(x = c(0.8, 1.2, 2), y = -5, paste("mouse", 1:3), xpd = T)
segments(x0 = c(0.7, 1.9), x1 = c(1.3, 2.1), y0 = -5.5, y1 = -5.5, xpd = T)
text(x = 1:2, y = -6, label = c("NY", "SY"), xpd = T)
nmouse = table(norm_per_coc$mouse[norm_per_coc$channel=="Esr2" & norm_per_coc$coc!="no"])
text(x = c(0.8,2, 1.2), y = -6.5, labels = paste0("n = ", nmouse), xpd= T)
segments(x0 = c(1, 0.8, 1, 2), x1 = c(2, 1.2, 1, 2), y0 = c(7.5, 7, 7, 7), y1 = c(7.5, 7, 7.5, 7.5), xpd = T)
text(x = 1.5, y = 7.7, label = paste("pval =", round(aov$`Pr(>Chisq)`[2],4)), xpd = T)
dev.off()

