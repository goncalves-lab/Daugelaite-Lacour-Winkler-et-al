#produces all files which are necessary for CNV analysis

library("infercnv")
library("Rgb")
library("stringr")

ann_db <- read.gtf("gencode.vM28.annotation.gtf")
ann_db <- ann_db[ann_db$feature=="gene",]
ann_db <- ann_db[,c(1,4,5,9,11)]
ann_db$gene_id <- str_replace(ann_db$gene_id,
               pattern = ".[0-9]+$",
               replacement = "")
ann_db <- ann_db[,c(4,1,2,3)]

load("gc_class.Rdata") #produced using classifier scripts
bad <- names(gc_class)[gc_class=="YS_bad"]
bad <- gsub("gc", "e", bad)

load("IVF.RData") #produced in scenic.R script

good <- names(gc_class)[gc_class=="YS_good"]
good <- gsub("gc", "e", good)

bad <- data.frame(label=bad, condition="bad")
good <- data.frame(label=good, condition="good")
metadata <- rbind(bad, good)

exprs_matrix <- IVF[["RNA"]]@counts

metadata_control <- metadata[metadata$condition=="good",]
metadata_so_bad <- metadata[metadata$condition=="bad",]


metadata_control$description <- "normal"
metadata_so_bad$description <- "malignant_1"

metadata <- rbind(metadata_control, metadata_so_bad)
exprs_matrix <- exprs_matrix[,colnames(exprs_matrix)%in%metadata$label]

metadata <- metadata[,c(1,3)]
metadata <- metadata[metadata$label%in%colnames(exprs_matrix),]
exprs_matrix <- exprs_matrix[,match(metadata$label, colnames(exprs_matrix))]

ann_db <- ann_db[ann_db$gene_id%in%rownames(exprs_matrix),]
exprs_matrix <- exprs_matrix[rownames(exprs_matrix)%in%ann_db$gene_id,]
exprs_matrix <- exprs_matrix[match(ann_db$gene_id,rownames(exprs_matrix)),]
write.table(as.data.frame(exprs_matrix), "exprs_matrix.tsv", quote=F, row.names = T, col.names = T,  sep = "\t")
write.table(metadata, "metadata.tsv", quote=F, row.names = F, col.names = F,  sep = "\t")
write.table(ann_db, "ann_db.tsv", quote=F, col.names = F, row.names = F, sep = "\t")


#######################################
# create the infercnv object

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="exprs_matrix.tsv",
                                    annotations_file="metadata.tsv",
                                    delim="\t",
                                    gene_order_file="ann_db.tsv",
                                    ref_group_names=c("normal"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir_sub",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T, HMM_report_by="cell",
                             HMM=T, analysis_mode='cells'
)

