
rm(list = ls())
memory.limit()
memory.limit(40000)
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)

AD.microglia<- merge(AD,y = c(ADcKO),add.cell.ids =NULL, project = "AD.microglia")
#Perform integration
ifnb.list <- SplitObject(AD.microglia, split.by = "orig.ident")
#length(ifnb.list)
#ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = ifnb.list)
AD.microglia.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
save(AD.microglia.anchors,file = 'AD.microglia.anchors.Rdata')

AD.microglia.anchors
anchorset = AD.microglia.anchors
# this command creates an 'integrated' data assay
AD.microglia.combined <- IntegrateData(anchorset = AD.microglia.anchors,dims = 1:30, k.weight = 90)


save(AD.microglia.combined,file = 'AD.microglia.combined.Rdata')

#Perform an integrated analysis 甯歌鍒嗘瀽鍏堜竴???
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(AD.microglia.combined) <- "integrated"

VlnPlot(AD.microglia.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4,group.by="orig.ident")
plot1 <- FeatureScatter(AD.microglia.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="orig.ident")
plot2 <- FeatureScatter(AD.microglia.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
plot1 + plot2 #榛樿group.by=cluster锛屽彲浠ヨ缃垚"groups"

# Run the standard workflow for visualization and clustering
AD.microglia.combined <- ScaleData(AD.microglia.combined,verbose = FALSE)
AD.microglia.combined <- RunPCA(AD.microglia.combined, npcs = 50, verbose = FALSE)
DimPlot(AD.microglia.combined, reduction = "pca")
ElbowPlot(AD.microglia.combined,ndims = 60)
#microglia.combined <- JackStraw(microglia.combined, num.replicate = 100)
#microglia.combined <- ScoreJackStraw(microglia.combined, dims = 1:20)
#JackStrawPlot(microglia.combined, dims = 1:20)

AD.microglia.combined <- RunUMAP(AD.microglia.combined, reduction = "pca", dims = 1:30)
AD.microglia.combined <- FindNeighbors(AD.microglia.combined, reduction = "pca", dims = 1:30)
AD.microglia.combined <- FindClusters(AD.microglia.combined, resolution = 1)

# Visualization
DimPlot(AD.microglia.combined, reduction = "umap",label =TRUE,
        label.size = 4,
        label.color = "black")
p1 <- DimPlot(AD.microglia.combined, reduction = "umap", split.by = "orig.ident",label =TRUE,
              label.size = 4,
              label.color = "black")
p2 <- DimPlot(AD.microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(AD.microglia.combined, reduction = "umap", group.by = "orig.ident") #microglia澶ц嚧鍒嗙兢
DimPlot(AD.microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)  #姣忕粍鍒嗙兢澶т綋
DimPlot(AD.microglia.combined, reduction = "umap",  label = TRUE,split.by = "orig.ident")

DefaultAssay(AD.microglia.combined) <- "RNA"
cluster2.markers <- FindMarkers(AD.microglia.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

AD.microglia.combined.markers <- FindAllMarkers(AD.microglia.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AD.microglia.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

AD.microglia.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(AD.microglia.combined) <- "integrated"
DoHeatmap(AD.microglia.combined, features = top10$gene) + NoLegend()


saveRDS(AD.microglia.combined.markers, file = "AD.microglia.combined_Markers.rds")

cluster0.markers <- FindMarkers(microglia.combined, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(microglia.combined, features = c("Crybb1", "Fscn1"))
VlnPlot(microglia.combined, features = c("Crybb1", "Fscn1"), slot = "counts", log = TRUE)
FeaturePlot(microglia.combined, features = c("Crybb1", "Fscn1"))

saveRDS(AD.microglia.combined, file = "AD.microglia.combined_Umap.rds")


#microglia提取后再分群
DefaultAssay(AD.microglia.combined_Umap) <- "RNA"
FeaturePlot(AD.microglia.combined_Umap, features = c("Tmem119","P2ry12","Trem2","Cx3cr1","Hexb"), min.cutoff = "q9",split.by = "orig.ident")
FeaturePlot(AD.microglia.combined_Umap, features = c("Ly6c2","Ly6c1"), min.cutoff = "q9",split.by = "orig.ident")
VlnPlot(AD.microglia.combined_Umap, features = c("Tmem119","P2ry12","Trem2","Cx3cr1","Hexb"), split.by = "seurat_clusters",
        pt.size = 0, combine = FALSE)
#VlnPlot(AD.microglia.combined_Umap, features = c("Plp1","Mbp","St18"), split.by = "seurat_clusters", pt.size = 0, combine = FALSE)

sub.AD.microglia <- subset(AD.microglia.combined_Umap, idents = c("0","1","3","4","5","7","10","12","13","14","17","18","15"))
sub.AD.microglia[["percent.rp"]]  = PercentageFeatureSet(sub.AD.microglia, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(sub.AD.microglia, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")

VlnPlot(sub.AD.microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
#sub.microglia <- subset(sub.microglia, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.rp < 15)

#ifnb.list <- SplitObject(sub.microglia, split.by = "orig.ident")
#length(ifnb.list)
#ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#x <- NormalizeData(x)
#x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
# select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = ifnb.list)
#sub_microglia.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
#save(sub_microglia.anchors,file = 'sub_microglia.anchors.Rdata')
#sub_microglia.anchors
#anchorset = sub_microglia.anchors
# this command creates an 'integrated' data assay
#sub.microglia <- IntegrateData(anchorset = sub_microglia.anchors,dims = 1:30, k.weight = 90)
#save(sub.microglia,file = 'sub.microglia.Rdata')

#Perform an integrated analysis 甯歌鍒嗘瀽鍏堜竴???
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay


VlnPlot(sub.AD.microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4,group.by="orig.ident")
#plot1 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="orig.ident")
#plot2 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
#plot1 + plot2 #榛樿group.by=cluster锛屽彲浠ヨ缃垚"groups"
sub.AD.microglia <- NormalizeData(sub.AD.microglia, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sub.AD.microglia<- FindVariableFeatures(sub.AD.microglia, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sub.AD.microglia), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sub.AD.microglia)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2


sub.AD.microglia <- ScaleData(sub.AD.microglia, verbose = FALSE)
sub.AD.microglia <- RunPCA(sub.AD.microglia, features = VariableFeatures(object = sub.AD.microglia))

#sub.microglia <- RunPCA(sub.microglia, npcs = 30, verbose = FALSE)
sub.AD.microglia <- RunUMAP(sub.AD.microglia, reduction = "pca", dims = 1:30)
sub.AD.microglia <- FindNeighbors(sub.AD.microglia, reduction = "pca", dims = 1:30)
sub.AD.microglia <- FindClusters(sub.AD.microglia, resolution = 1.4)

DimPlot(sub.AD.microglia, reduction = "umap",label = TRUE)
p1 <- DimPlot(sub.AD.microglia, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sub.AD.microglia, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(sub.AD.microglia, reduction = "umap", split.by = "orig.ident",label = TRUE)

DefaultAssay(sub.AD.microglia) <- "RNA"
sub.AD.microglia.markers <- FindAllMarkers(sub.AD.microglia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sub.AD.microglia.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

sub.AD.microglia.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(sub.AD.microglia,features = top10$gene) + NoLegend()

saveRDS(sub.AD.microglia,file = 'sub.AD.microglia_Umap.rds')
saveRDS(sub.AD.microglia.markers,file = 'sub.microglia.markers.rds')

AD.cluster9.markers <- FindMarkers(sub.AD.microglia_Umap, ident.1 = 9, min.pct = 0.25)
#head(cluster10.markers, n = 5)


#查看细胞分群来源比例
library(gplots)
tab.1=table(sub.AD.microglia_Umap$orig.ident,sub.AD.microglia_Umap$seurat_clusters) 
balloonplot(tab.1)

FeaturePlot(sub.microglia_Umap, features = c("Tmem119","P2ry12","Trem2"), min.cutoff = "q9",split.by = "orig.ident")
VlnPlot(sub.AD.microglia_Umap, features = c("Folr1"), split.by = "seurat_clusters",group.by = "orig.ident",
        pt.size = 0, combine = FALSE)



#以下代码未跑
#以下代码未跑
#检测细胞周期分群影响
#microglia.combined_Umap.rds
library(Seurat)
library(tidyverse)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
microglia.combined_Umap <- CellCycleScoring(microglia.combined_Umap, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(microglia.combined_Umap@meta.data,2)
DimPlot(microglia.combined_Umap,reduction = "umap")

microglia.combined_Umap <- SetIdent(microglia.combined_Umap, value = "seurat_clusters")
saveRDS(microglia.combined_Umap, file = "microglia.combined_Umap.rds")













