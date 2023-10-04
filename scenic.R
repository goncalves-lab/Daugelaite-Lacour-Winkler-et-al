#performs SCENIC analysis
library(SCENIC)
library(gsubfn) 
library(superheat)
library(RColorBrewer)
library(biomaRt)

#setting SCENIC setting and filtering expr. matrix
SCENIC_organizeR <- function(dbDir, myDatasetTitle, exprs_matrix){
  org <- "mgi" # or hgnc, or dmel
  dbDir <- dbDir # RcisTarget databases location
  myDatasetTitle <- myDatasetTitle # choose a name for your analysis
  data(defaultDbNames)
  dbs <- defaultDbNames[[org]]
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=6)
# (Adjust minimum values according to your dataset)
  genesKept <- geneFiltering(exprs_matrix, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprs_matrix),
                           minSamples=ncol(exprs_matrix)*.01)
  exprMat_filtered <- exprs_matrix[genesKept, ]
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered <- log2(exprMat_filtered+1)
  list(exprMat_filtered, scenicOptions)}

#running Genie3 (Building the gene regulatory network (GRN))
SCENIC_runner <- function(exprMat_filtered, scenicOptions){
  set.seed(123)
  runGenie3(exprMat_filtered, scenicOptions)
  scenicOptions@settings$seed <- 123
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

#Select potential direct-binding targets (regulons) based on DNA-motif analysis (RcisTarget: TF motif analysis)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

#Analyzing the network activity in each individual cell (AUCell)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  cell_info <- cell_info[cell_info$label%in%colnames(getAUC(regulonAUC)),]
  regulonActivity_byCellType <- sapply(split(cell_info$label, cell_info$cells),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  list(regulonAUC, regulonActivity_byCellType)}

SCENIC_binarizeR <- function(scenicOptions, cell_info){
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  minPerc <- .7
  binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
  cellInfo_binarizedCells <- cell_info[which(rownames(cell_info)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
  regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$label), 
                                                 function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
  list(binaryRegulonActivity, regulonActivity_byCellType_Binarized, cellInfo_binarizedCells)
}

ploteR <- function(auc_matrix, col_vec){
  if (length(col_vec)>2){
    heat <- c(0, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8,1)
    #breaks=c(0, 0.2, 0.4, 0.6,0.8)
  } 
  else {
    heat <- c(0, 1)
    breaks <- c(0,1)
  }
  print(superheat(auc_matrix,  
            #grid.hline.col = "white",
            #grid.vline.col = "white",
            
            # rotate bottom label text
            bottom.label.text.angle = 90,
            heat.pal = col_vec, heat.pal.values = heat,
            #legend.breaks = breaks, smooth.heat = F,
            grid.hline.size = 1, grid.vline.size=1,
            scale=F, left.label.col = "white", bottom.label.col = "white",left.label.text.size = 7, bottom.label.size = 0.5, pretty.order.cols = T,
            pretty.order.rows = T, col.dendrogram = T,  force.bottom.label = TRUE, dist.method="euclidean"))
}
########################################################

#loading expression matrix (wxprs_matrix) and metadata (cell_info)

exprs_matrix <- read.csv("reads_all_samples.csv", row.names = 1)
colnames(exprs_matrix) <- gsub("_rep", "", colnames(exprs_matrix))
#exprs_matrix <- t(exprs_matrix)
cell_info <- data.frame(label = colnames(exprs_matrix))
cell_info$label <- gsub("_rep", "", cell_info$label)
cell_info$cells <- substr(cell_info$label,nchar(cell_info$label)-3, nchar(cell_info$label)-2)
cell_info$day <- substr(cell_info$label,1, nchar(cell_info$label)-4)
cell_info$nGene <- colSums(exprs_matrix >0)

exprs_matrix <- make.names(exprs_matrix, unique = T)
rownames(cell_info) <- make.names(cell_info$label, unique = T)
ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensembl_list=getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values = row.names(exprs_matrix), mart=ensembl)

for(i in seq(from=1, to=nrow(ensembl_list))){
  if (ensembl_list[i,2]==""){
    ensembl_list[i,2] = "no_id"
    }
}
full_list=merge(exprs_matrix, ensembl_list,by.x=0,by.y="ensembl_gene_id")
rownames(full_list) <- make.names(full_list$external_gene_name, unique=T)
exprs_matrix <- as.matrix(full_list[,2:221])

list[exprMat_filtered, scenicOptions] <- SCENIC_organizeR("cis_Target_db", "SCENIC OCs & GCs", exprs_matrix)
list[regulonAUC, regulonActivity_byCellType] <- SCENIC_runner(exprMat_filtered, scenicOptions)

aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
list[binaryRegulonActivity, regulonActivity_byCellType_Binarized, cellInfo_binarizedCells] <- SCENIC_binarizeR(scenicOptions, cell_info)


auc_matrix <- getAUC(regulonAUC)
auc_matrix <- auc_matrix[rowSums(auc_matrix)> 0.2,] 
auc_matrix_oc <- auc_matrix[, grep(".*OC.*", colnames(auc_matrix))]
auc_matrix_gc <- auc_matrix[, grep(".*GC.*", colnames(auc_matrix))]

write.csv(auc_matrix_oc, "auc_matrix_oc", quote = F)
write.csv(auc_matrix_gc, "auc_matrix_gc", quote = F)

