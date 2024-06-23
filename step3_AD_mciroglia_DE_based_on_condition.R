
library(Seurat)
library(cowplot)

DefaultAssay(sub.microglia_Umap) <- "RNA"

obj.list <- SplitObject(sub.AD.microglia_Umap, split.by = "orig.ident")
AD.microglia<-obj.list$AD
ADcKO.microglia<-obj.list$ADcKO
AD.microglia$stim <- "AD"
ADcKO.microglia$stim <- "ADcKO"

#SplitObject(sub.AD.microglia_Umap, split.by = "orig.ident")
Idents(sub.AD.microglia_Umap) <- "orig.ident"
ADvsADcKO <- FindMarkers(sub.AD.microglia_Umap, ident.1 = "AD", ident.2 = "ADcKO", verbose = FALSE)

###
FeaturePlot(sub.AD.microglia_Umap, features = c("Ttr"), split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(sub.AD.microglia_Umap, features = c("Ttr"), split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

FeaturePlot(sub.microglia_Umap, features = c("Nr3c2","Nr3c1"), split.by = "orig.ident", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(sub.microglia_Umap, features = c("Nr3c2","Nr3c1"), split.by = "orig.ident", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

plots <- VlnPlot(sub.microglia_Umap, features = c("Fkbp5","Ttr"), split.by = "orig.ident",pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
wrap_plots(plots = plots, ncol = 1)

a<-as.data.frame(ADvsADcKO, genes = Seurat::VariableFeatures(x), fix_names = TRUE)
write.csv(a,file="sub.AD.microglia_ADvsADcKO_DE_genes.csv")





