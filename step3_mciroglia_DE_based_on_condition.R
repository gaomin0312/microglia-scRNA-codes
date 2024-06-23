
library(Seurat)
library(cowplot)

DefaultAssay(sub.microglia_Umap) <- "RNA"

obj.list <- SplitObject(sub.AD.microglia_012, split.by = "orig.ident")
cKO.microglia<-obj.list$cKO
Flox.microglia<-obj.list$Flox
Flox.microglia$stim <- "ctrl"
cKO.microglia$stim <- "cKO"

#SplitObject(sub.microglia_Umap, split.by = "orig.ident")
Idents(sub.AD.microglia_Umap) <- "orig.ident"
sub.AD.microglia_012 <- subset(AD.microglia.combined_Umap, idents = c("0","1","2"))
ADvsADcKO_012 <- FindMarkers(sub.AD.microglia_012, ident.1 = "AD", ident.2 = "ADcKO", verbose = FALSE)

FeaturePlot(sub.microglia_Umap, features = c("Cst7"), split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(sub.microglia_Umap, features = c("Tmem119"), split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
#wrap_plots(plots = plots, ncol = 1)
plots <- VlnPlot(sub.microglia_Umap, features = c("H2-D1"),pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

a<-as.data.frame(cKOvsCtrl, genes = Seurat::VariableFeatures(x), fix_names = TRUE)
write.csv(a,file="sub.microglia_cKOvsFlox_DE_genes.csv")

cKO_cluster9 <- subset(sub.microglia_Umap, idents = c("9"))
FeaturePlot(cKO_cluster9, features = c("P2ry12","Tmem119","Mrc1","Lyz2","Cst7"), split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(sub.microglia_Umap, features = c("Tmem119"), split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)


rownames(avg.sub.microglia )




