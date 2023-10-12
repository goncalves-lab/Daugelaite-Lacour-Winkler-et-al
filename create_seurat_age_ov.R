#Generating the Seurat object from the count table and metadata file

library(Seurat)

#Import files from ArrayExpress
reads_df = read.csv("count_table_all_samples.csv", header = T, row.names = 1)
metafilter = read.csv("metadata_filtered.csv") #sample metadata table


sample_vec =  metafilter$Name
reads_df = reads_df[,sample_vec]

#Create Seurat object
ovulation = CreateSeuratObject(counts = reads_df, project = "ovulation")
#Running normalisation and dimensionality reductions
ovulation = SCTransform(ovulation, verbose = FALSE,return.only.var.genes = FALSE)
ovulation = RunPCA(ovulation, features = VariableFeatures(object = ovulation))
ovulation = RunUMAP(ovulation, dims = 1:30)
ovulation = FindNeighbors(ovulation, dims = 1:30)
ovulation = FindClusters(ovulation, resolution = 0.5)

label_df = as.data.frame(names(ovulation$orig.ident))
label_df$label = unlist(lapply(strsplit(label_df$`names(ovulation$orig.ident)`, "C"), function(X){paste0(X[1],"C")}))

#Parsing the sample names to get age, ovulation and cell type
list_plot = list()
for (i in seq(from=1, to=nrow(label_df))){
  lab_vector = label_df[i,2]
  if (grepl("S[0O]",label_df[i,2])){
    lab_vector[2] = "SO"
  } else {lab_vector[2] = "NO"}
  if (grepl("MO",label_df[i,2],fixed=T)){
    lab_vector[3] = "OC"
  } else {lab_vector[3] = "GC"}
  if (grepl("3M",label_df[i,2],fixed=T)){
    lab_vector[4] = "Y"
  } else {lab_vector[4] = "O"}
  list_plot[[i]] = lab_vector
}

label_df_ex = data.frame(matrix(unlist(list_plot), nrow=length(list_plot), byrow=T))
colnames(label_df_ex) = c("label", "ovulation", "cell_type", "age")
label_df_ex$mouse = metafilter$Individual[match(names(ovulation$orig.ident), metafilter$Name)]

#Adding age, ovulation and cell type to the Seurat object
ovulation$label = label_df_ex$label
ovulation$age = label_df_ex$age
ovulation$cell = label_df_ex$cell_type
ovulation$ovulation = label_df_ex$ovulation
ovulation$mouse = label_df_ex$mouse

save(ovulation, file="ovulation.Rdata")
