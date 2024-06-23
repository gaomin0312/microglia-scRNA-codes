rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

library(dplyr)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
getwd()
setwd("D:/0scRNA_data_analysis/microglia_cKO_Flox/cKO8mo-male-FACS_microglia/outs/filtered_feature_bc_matrix-bySeurat") #gm 设置目录需要带�?"D:/..."
# ############################################################################################################
# ## load data 
# ############################################################################################################


dir.10x = 'D:/0scRNA_data_analysis/0-data-cKO8mo-male-MI-filtered_feature_bc_matrix-bySeurat/'
list.files(dir.10x)

cKO_MI.data <- Read10X(data.dir=dir.10x)
# Set up
cKO <- CreateSeuratObject(counts = cKO_MI.data, project = "cKO",
                          min.cells = 3, min.features = 200) #过滤参数

cKO[["percent.mt"]] <- PercentageFeatureSet(cKO, pattern = "^mt-")
cKO[["percent.rp"]]  = PercentageFeatureSet(cKO, pattern = "^Rp[sl][[:digit:]]")

head(cKO@meta.data, 5)
VlnPlot(cKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cKO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cKO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cKO <- subset(cKO, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
cKO;#AD;WT

setwd("D:/0scRNA_data_analysis/microglia_cKO_Flox/FloxCtrl8mo-male-FACS_microglia/outs/filtered_feature_bc_matrix") #gm 设置目录需要带�?"D:/..."
dir.10x = 'D:/0scRNA_data_analysis/microglia_cKO_Flox/FloxCtrl8mo-male-FACS_microglia/outs/filtered_feature_bc_matrix/'
list.files(dir.10x)

Flox_MI.data <- Read10X(data.dir=dir.10x)
# Set up
Flox <- CreateSeuratObject(counts =Flox_MI.data, project = "Flox",
                          min.cells = 3, min.features = 200) #过滤参数

Flox[["percent.mt"]] <- PercentageFeatureSet(Flox, pattern = "^mt-")
Flox[["percent.rp"]]  = PercentageFeatureSet(Flox, pattern = "^Rp[sl][[:digit:]]")
head(Flox@meta.data, 5)
VlnPlot(Flox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Flox, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Flox, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Flox <- subset(Flox, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
Flox;#AD;WT

setwd("C:/Users/gaomi/Desktop/microglia_AD_ADcKO_cKO_Flox/AD8m/outs") #gm 设置目录需要带�?"D:/..."
# ############################################################################################################
# ## load data 
# ############################################################################################################


dir.10x = 'C:/Users/gaomi/Desktop/microglia_AD_ADcKO_cKO_Flox/AD8m/outs/'
list.files(dir.10x)

AD_MI.data <- Read10X(data.dir=dir.10x)
# Set up
AD <- CreateSeuratObject(counts = AD_MI.data, project = "AD",
                          min.cells = 3, min.features = 200) #过滤参数

AD[["percent.mt"]] <- PercentageFeatureSet(AD, pattern = "^mt-")
AD[["percent.rp"]]  = PercentageFeatureSet(AD, pattern = "^Rp[sl][[:digit:]]")

head(AD@meta.data, 5)
VlnPlot(AD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(AD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
AD <- subset(AD, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
AD;#AD;WT

setwd("C:/Users/gaomi/Desktop/microglia_AD_ADcKO_cKO_Flox/ADcKO8m/outs") #gm 设置目录需要带�?"D:/..."
# ############################################################################################################
# ## load data 
# ############################################################################################################


dir.10x = 'C:/Users/gaomi/Desktop/microglia_AD_ADcKO_cKO_Flox/ADcKO8m/outs/'
list.files(dir.10x)

ADcKO_MI.data <- Read10X(data.dir=dir.10x)
# Set up
ADcKO <- CreateSeuratObject(counts = ADcKO_MI.data, project = "ADcKO",
                         min.cells = 3, min.features = 200) #过滤参数

ADcKO[["percent.mt"]] <- PercentageFeatureSet(ADcKO, pattern = "^mt-")
ADcKO[["percent.rp"]]  = PercentageFeatureSet(ADcKO, pattern = "^Rp[sl][[:digit:]]")

head(ADcKO@meta.data, 5)
VlnPlot(ADcKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ADcKO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ADcKO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ADcKO <- subset(ADcKO, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
ADcKO;#AD;WT

AD[["percent.rp"]]  = PercentageFeatureSet(AD, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(AD, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")
VlnPlot(AD, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
AD <- subset(AD, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA<20000 & percent.rp < 25 & percent.mt < 5)

ADcKO[["percent.rp"]]  = PercentageFeatureSet(ADcKO, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(ADcKO, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")
VlnPlot(ADcKO, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
ADcKO <- subset(ADcKO, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA<20000 & percent.rp < 25 & percent.mt < 5)

cKO[["percent.rp"]]  = PercentageFeatureSet(cKO, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(cKO, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")
VlnPlot(cKO, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
cKO <- subset(cKO, subset = nFeature_RNA > 200 & nFeature_RNA < 4500  & nCount_RNA<15000 & percent.rp < 20 & percent.mt < 5)

Flox[["percent.rp"]]  = PercentageFeatureSet(Flox, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(Flox, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")
VlnPlot(Flox, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
Flox <- subset(Flox, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA<15000 & percent.rp < 20 & percent.mt < 5)


# Set up
cKO <- NormalizeData(cKO, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#cKO <- NormalizeData(cKO)
cKO <- FindVariableFeatures(cKO, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Flox <- NormalizeData(Flox, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#Flox <- NormalizeData(Flox)
Flox <- FindVariableFeatures(Flox, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

AD <- NormalizeData(AD, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#cKO <- NormalizeData(cKO)
AD <- FindVariableFeatures(AD, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

ADcKO <- NormalizeData(ADcKO, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#cKO <- NormalizeData(cKO)
ADcKO <- FindVariableFeatures(ADcKO, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top10_1 <- head(VariableFeatures(cKO), 10)
top10_2 <- head(VariableFeatures(Flox), 10)
top10_3 <- head(VariableFeatures(AD), 10)
top10_4 <- head(VariableFeatures(ADcKO), 10)

setwd("C:/Users/gaomi/Desktop/microglia_AD_ADcKO_cKO_Flox")
saveRDS(cKO, file = "cKO.rds")
saveRDS(Flox, file = "Flox.rds")
saveRDS(AD, file = "AD.rds")
saveRDS(ADcKO, file = "ADcKO.rds")
