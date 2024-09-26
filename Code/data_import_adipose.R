library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)
library(scRepertoire)

## set working directory
setwd("/Users/chengyuye/Library/CloudStorage/OneDrive-EmoryUniversity/Rotation_2_Li Lab/scRNA_TCR/Data/p22285-s001_Adipose-Tregs_WT-KO-mix")

### read in data
adipo <- Read10X(data.dir = './count/sample_filtered_feature_bc_matrix/')
adipo_hto <- adipo$`Antibody Capture` %>% as.data.frame()

## check out samples
rownames(adipo_hto)

## remove the "Hash3Mouse_C0303_TotalSeqC" 
adipo_hto <- adipo_hto[-3, ]

table(colSums(adipo_hto)== 0)
adipo_hto <- adipo_hto[colSums(adipo_hto)>0]
dim(adipo_hto)


## check out the gene expression matrix
adipo_gex <- adipo$`Gene Expression`

#  get the intersection
joint.bcs <- intersect(colnames(adipo_gex), colnames(adipo_hto))
length(joint.bcs)

# Subset RNA and HTO counts by joint cell barcodes
ear_gex <- ear_gex[, joint.bcs]
ear_hto <- as.matrix(ear_hto[, joint.bcs])

# Setup Seurat object
adipo.hashtag <- CreateSeuratObject(counts = adipo_gex)

# Normalize RNA data with log normalization
adipo.hashtag <- NormalizeData(adipo.hashtag)
# Find and scale variable features
adipo.hashtag <- FindVariableFeatures(adipo.hashtag, selection.method = "mean.var.plot")
adipo.hashtag <- ScaleData(adipo.hashtag, features = VariableFeatures(adipo.hashtag))
adipo.hashtag

# Add HTO data as a new assay independent from RNA
adipo.hashtag[["HTO"]] <- CreateAssayObject(counts = adipo_hto)
names(adipo.hashtag)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
adipo.hashtag <- NormalizeData(adipo.hashtag, 
                               assay = "HTO", 
                               normalization.method = "CLR")


head(adipo.hashtag@meta.data)

VlnPlot(adipo.hashtag, features = c("nFeature_RNA", "nCount_RNA",
                                   "nCount_HTO","nFeature_HTO"), 
        ncol = 2,pt.size = 0)


## demultiplexing
adipo.hashtag

adipo.hashtag <- HTODemux(adipo.hashtag,
                          assay = "HTO", 
                          positive.quantile = 0.99)
adipo.hashtag
table(adipo.hashtag$HTO_classification.global)
adipo.hashtag@assays
table(adipo.hashtag$hash.ID)

# Doublet Negative  Singlet 
# 504      654     7878

## remove nagative cells
adipo.hashtag <- subset(adipo.hashtag, 
                        subset = HTO_classification.global != 'Negative')


## split the sample based on Ab
### WT ---- 4213 cells before QC
adipo_Hash1_C0301 <- subset(adipo.hashtag,
                            subset = HTO_maxID == 'Hash1Mouse-C0301-TotalSeqC')


### Srebf2KO ---- 4169 cells before QC
adipo_Hash2_C0302 <- subset(adipo.hashtag,
                            subset = HTO_maxID == 'Hash2Mouse-C0302-TotalSeqC')

## add metadata
adipo_Hash1_C0301$Group <- 'WT'

adipo_Hash2_C0302$Group <- 'Srebf2KO'



### save the data
saveRDS(adipo.hashtag, '../processed_data/combined_hashtag_negative_removed.rds')
saveRDS(adipo_Hash1_C0301, '../processed_data/WT_adipose.rds')
saveRDS(adipo_Hash2_C0302, '../processed_data/Srebf2KO_adipose.rds')


