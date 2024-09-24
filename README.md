## Citation

If you use our work, please cite: 
```
Daugelaite K, Lacour P, Winkler I, Koch M, Schneider A, Schneider N, Tolkachov A, Nguyen XP, Vilkaite A, Rehnitz J, Odom DT, Goncalves A. (2023) Superovulation and ageing perturb oocyte-granulosa cell transcriptomes and communication. bioRxiv (2023). doi: https://doi.org/10.1101/2023.10.30.563978
```

## From raw counts to R objects


`create_seurat_age_ov.R` - for natural and superovulated, young and old oocytes and granulosa cells Smart-Seq2 data (E-MTAB-13479)  
`create_seurat_totalrna.R` - for natural and superovulated oocytes total-RNA seq data (E-MTAB-13474)  
`create_seurat_ivf_mouse.R` - for IVF-derived mouse embryos (morula or blastocyst) and corresponding granulosa cells, Smart-seq2 data (E-MTAB-13480)

These scripts create the Seurat objects used by the other scripts from the raw count tables.

`Scnorm.R` - normalizes count data using the SCnorm method to take into account gene length (used for cell communication and classifier scripts).


## Differential gene expression and overrepresentation

`dge.R` - differential expression analysis using DESeq2 for aging and superovulation dataset  
`ora.R` - over-representation analysis of genes found by DESeq2

`dge_SNvS.R` - differential expression analysis using DESeq2 between S and SN granulosa cells 
(as identified by cell-to-cell communication analysis and transcription factor activities)



## Total-RNA analysis (repolyadenylation, deadenylation, and degradation)

`totalrna_vs_smartseq.R` - compares the expression of known genes between natural and superovulated oocytes
in a polyA-biased technology (Smart-Seq2) and a non-biased one (total RNA)

## Cell-to-cell communication

`cell_communication.R` - computes ligand-receptor interaction score based on gene expression level and CellChatDB annotation

## SCENIC and AUCell

`scenic.R` - runs SCENIC analysis on oocytes and granulosa cells from the aging and superovulation dataset  
`scenic_post.R` - tests for significant differentially active pathways between conditions  
`tf_scenic_pathway.R` - computes the overlap between the TFs targets and the pathways, plots the results in a heatmap

`aucell.R` - computes pathway activity scores

## Granulosa cells classifier

`data_preparation_genes.R` - selects genes that will be used in the gene classifier 
(based on differentially expressed genes (DEG) between S and SN granulosa cells)  
`data_preparation_tfs.R` - selects genes that will be used in the TF classifier (based on SCENIC results)  
`auc_classifier.R` - trains different granulosa classifiers using TF activity scores  
`genes_classifier.R` - trains different granulosa classifiers using DEG  
`gc_scenic_scoring_classifier.R` - computes TF activity scores of new samples using the same regulons as the ones in the training dataset 
(results from the SCENIC analysis)  
`classifier_combined.R` - predicts the class of new granulosa cells using the two classifiers

## Embryo development and copy number variation

`pseudotime_embryos.R` - creates a reference developmental trajectory and calculates a developmental pseudotime for each embryo
to assess link between granulosa cell classification and developmental transcriptional trajectory
`CNV_prep.R` - prepares embryo data for inferCNV run  
`CNV_runner.R` - runs inferCNV on embryo data

## Validation experiments (HCR and qPCR)

`hcr_analysis_and_plots.R` - validation of Esr2 expression in natural and superovulated young granulosa cells using HCR fluorescence  
`qPCR_analysis_and_plots.R` - qPCR quantification of genes used in the granulosa cells classifier

## Shannon entropy

`Shannon_entropy.R` - computes differential Shannon entropy for the aging and superovulation dataset

## Human data analysis

`human_dge_gsea.R` - computes differential gene expression on human granulosa cells (E-MTAB-13496) and 
compares the enriched pathways identified using fgsea to the ones found in mouse

## Interaction between aging and superovulation

`pca_projection.R` - uses a PCA projection approach to summarize the non-linearity between aging and superovulation effects

## Old scripts used in previous versions of the manuscript

`pseudotime_oc_gc.R` - performs pseudotime analysis based on highly variable genes or pathways of interest (e.g. meiosis)

