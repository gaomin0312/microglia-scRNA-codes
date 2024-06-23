
rm(list = ls())
memory.limit()
memory.limit(40000)
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)

Flox_1<-subset(Flox, downsample=4000)
microglia<- merge(cKO,y = c(Flox_1, AD, ADcKO),add.cell.ids =NULL, project = "microglia")
#Perform integration
ifnb.list <- SplitObject(microglia, split.by = "orig.ident")
#length(ifnb.list)
#ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = ifnb.list)
microglia.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
save(microglia.anchors,file = 'microglia.anchors.Rdata')

microglia.anchors
anchorset = microglia.anchors
# this command creates an 'integrated' data assay
microglia.combined <- IntegrateData(anchorset = microglia.anchors,dims = 1:30, k.weight = 90)


save(microglia.combined,file = 'microglia.combined.Rdata')

#Perform an integrated analysis 常规分析先一???
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(microglia.combined) <- "integrated"

VlnPlot(microglia.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4,group.by="orig.ident")
plot1 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="orig.ident")
plot2 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
plot1 + plot2 #默认group.by=cluster，可以设置成"groups"

# Run the standard workflow for visualization and clustering
microglia.combined <- ScaleData(microglia.combined,verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 50, verbose = FALSE)
DimPlot(microglia.combined, reduction = "pca")
ElbowPlot(microglia.combined,ndims = 60)
#microglia.combined <- JackStraw(microglia.combined, num.replicate = 100)
#microglia.combined <- ScoreJackStraw(microglia.combined, dims = 1:20)
#JackStrawPlot(microglia.combined, dims = 1:20)

microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:30)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:30)
microglia.combined <- FindClusters(microglia.combined, resolution = 1)

# Visualization
DimPlot(microglia.combined, reduction = "umap",label =TRUE,
        label.size = 4,
        label.color = "black")
p1 <- DimPlot(microglia.combined, reduction = "umap", split.by = "orig.ident",label =TRUE,
              label.size = 4,
              label.color = "black")
p2 <- DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(microglia.combined, reduction = "umap", group.by = "orig.ident") #microglia大致分群
DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)  #每组分群大体
DimPlot(microglia.combined, reduction = "umap",  label = TRUE,split.by = "orig.ident")

DefaultAssay(microglia.combined) <- "RNA"
cluster2.markers <- FindMarkers(microglia.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

microglia.combined.markers <- FindAllMarkers(microglia.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(microglia.combined.markers,file = 'microglia.combined.markers.csv')
microglia.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

microglia.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(microglia.combined) <- "integrated"
DoHeatmap(microglia.combined, features = top10$gene) + NoLegend()


saveRDS(microglia.combined.markers, file = "microglia.combined_Markers.rds")

cluster0.markers <- FindMarkers(microglia.combined, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(microglia.combined, features = c("Crybb1", "Fscn1"))
VlnPlot(microglia.combined, features = c("Crybb1", "Fscn1"), slot = "counts", log = TRUE)
FeaturePlot(microglia.combined, features = c("Crybb1", "Fscn1"))

saveRDS(microglia.combined_Umap, file = "microglia.combined_Umap.rds")


#microglia??ȡ???ٷ?Ⱥ
DefaultAssay(microglia.combined_Umap) <- "RNA"
FeaturePlot(microglia.combined_Umap, features = c("Tmem119","P2ry12","Trem2","Cx3cr1","Hexb"), min.cutoff = "q9",split.by = "orig.ident")
FeaturePlot(microglia.combined_Umap, features = c("Ly6c2","Ly6c1"), min.cutoff = "q9",split.by = "orig.ident")
VlnPlot(microglia.combined_Umap, features = c("Camk2a"), split.by = "seurat_clusters",
        pt.size = 0, combine = FALSE)
VlnPlot(microglia.combined_Umap, features = c("Plp1","Mbp","St18"), split.by = "seurat_clusters",
        pt.size = 0, combine = FALSE)

sub.microglia <- subset(microglia.combined_Umap, idents = c("0","1","2","3","4","5","7","10","14","15","16","17","18","20"))
sub.microglia[["percent.rp"]]  = PercentageFeatureSet(sub.microglia, pattern = "^Rp[sl][[:digit:]]")
FeatureScatter(sub.microglia, feature1 = "nCount_RNA", feature2 = "percent.rp", group.by ="orig.ident")

VlnPlot(sub.microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
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

#Perform an integrated analysis 常规分析先一???
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay


VlnPlot(sub.microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4,group.by="orig.ident")
#plot1 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="orig.ident")
#plot2 <- FeatureScatter(microglia.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
#plot1 + plot2 #默认group.by=cluster，可以设置成"groups"
sub.microglia <- NormalizeData(sub.microglia, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
sub.microglia<- FindVariableFeatures(sub.microglia, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sub.microglia), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sub.microglia)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2


sub.microglia <- ScaleData(sub.microglia, verbose = FALSE)
sub.microglia <- RunPCA(sub.microglia, features = VariableFeatures(object = sub.microglia))

#sub.microglia <- RunPCA(sub.microglia, npcs = 30, verbose = FALSE)
sub.microglia <- RunUMAP(sub.microglia, reduction = "pca", dims = 1:30)
sub.microglia <- FindNeighbors(sub.microglia, reduction = "pca", dims = 1:30)
sub.microglia_Umap <- FindClusters(sub.microglia_Umap, resolution = 0.8)

DimPlot(sub.microglia_Umap, reduction = "umap",label = TRUE)
p1 <- DimPlot(sub.microglia_Umap, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sub.microglia_Umap, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(sub.microglia_Umap, reduction = "umap", split.by = "orig.ident",label = TRUE)

DefaultAssay(sub.microglia_Umap) <- "RNA"
sub.microglia.markers <- FindAllMarkers(sub.microglia_Umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sub.microglia.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

sub.microglia.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(sub.microglia,features = top10$gene) + NoLegend()

saveRDS(sub.microglia_Umap,file = 'sub.microglia_Umap.rds')
saveRDS(sub.microglia.markers,file = 'sub.microglia.markers.rds')
write.csv(sub.microglia.markers,file = 'sub.microglia.markers.csv')
cluster14.markers <- FindMarkers(sub.microglia_Umap, ident.1 = 9, min.pct = 0.25)
#head(cluster10.markers, n = 5)


#?鿴ϸ????Ⱥ��Դ????
library(gplots)
tab.1=table(sub.microglia_Umap$orig.ident,sub.microglia_Umap$seurat_clusters) 
balloonplot(tab.1)

FeaturePlot(sub.microglia_Umap, features = c("Apoe"), min.cutoff = "q9",split.by = "orig.ident")
VlnPlot(sub.microglia_Umap, features = c("Apoe"), split.by = "seurat_clusters",group.by = "orig.ident",
        pt.size = 0, combine = FALSE)



#???



´???δ??
#???´???δ??
#????ϸ?????ڷ?ȺӰ??
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













