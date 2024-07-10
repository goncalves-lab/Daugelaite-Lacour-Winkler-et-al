#runs CNV analysis using InferCNV

library("infercnv")
#library("Rgb")
library("stringr")


infercnv_obj = CreateInfercnvObject(raw_counts_matrix="exprs_matrix.tsv",
                                    annotations_file="metadata.tsv",
                                    delim="\t",
                                    gene_order_file="ann_db.tsv",
                                    ref_group_names=c("natural_morula", "natural_blastocyst"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir_sub",  # dir is auto-created for storing outputs
                             denoise=T, HMM_report_by="cell",
                             HMM=T, analysis_mode='cells'
)

