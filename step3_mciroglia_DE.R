library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(data.table)

# For performing differential expression after integration, we switch back to the original
library(multtest)
library(metap)
# data
DefaultAssay(microglia.combined_Umap) <- "RNA"

#Microglia
FeaturePlot(microglia.combined_Umap, features = c("P2ry12" ,"Tmem119","Cx3cr1","Hexb"), min.cutoff = "q9",split.by = "orig.ident")
microglia.combined.markers <- FindAllMarkers(microglia.combined_Umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
microglia.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

FeaturePlot(microglia.combined_Umap, features = c("P2ry12" ,"Tmem119","Cx3cr1","Hexb"), min.cutoff = "q9",split.by = "seurat_clusters")
plots <- VlnPlot(microglia.combined_Umap, features = c("P2ry12" ,"Tmem119","Cx3cr1","Hexb"), split.by = "seurat_clusters", group.by = "orig.ident",
                 pt.size = 0, combine = FALSE)
library(ggplot2)
wrap_plots(plots = plots, ncol = 1)
#画图配色
library(paletteer)
col<-paletteer_c("scico::berlin", 10)[c(1,3,4,9,5,2,6,8,10)]
DotPlot(microglia.combined_Umap,  features = c("Ccl3","Ccl4","Ccl12","Tnfaip3","Myc"), cols = col, dot.scale = 8, split.by = "seurat_clusters") +
  RotatedAxis()
plots <- VlnPlot(microglia.combined_Umap, features = c("Ccl3","Ccl4","Ccl12","Tnfaip3","Myc"), split.by = "seurat_clusters", group.by = "orig.ident",                pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#test甯歌FindAllMarkers()
microglia.combined.markers <- FindAllMarkers(microglia.combined_Umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
microglia.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

microglia.combined.markers_noRibo<-microglia.combined.markers[!grepl("^Rp[sl]",microglia.combined.markers$gene, ignore.case = F),]

a<-microglia.combined.markers_noRibo[order(microglia.combined.markers_noRibo$avg_log2FC),]
write.table(a,file = "microglia_markers_for_GSEA.xls")

top20<-microglia.combined.markers_noRibo%>%group_by(cluster)%>%top_n(n=20,wt=avg_log2FC)
#microglia.combined_Umap@assays$RNA@scale.data <- scale(microglia.combined_Umap@assays$RNA@data, scale = TRUE)
DoHeatmap(microglia.combined_Umap, features = top20$gene,size = 0.5,slot = "scale.data") + NoLegend()
#DoHeatmap(microglia.combined_Umap, features = rownames(microglia.combined.markers)) 
saveRDS(microglia.combined.markers_noRibo, file = "microglia.combined.markers_noRibo.rds")



#差异基因火山图
##cluster1vs0
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
deg0vs1<-FindMarkers(microglia.combined_Umap,ident.1 = '0',
                  ident.2 = '1')
deg0vs1_noRibo<-deg0vs1[!grepl("^Rp[sl]",rownames(deg0vs1), ignore.case = F),]

deg0vs1_noRibo[order(deg0vs1_noRibo$p_val),]
library(EnhancedVolcano) 
EnhancedVolcano(deg0vs1_noRibo,ylim = c(0, 400),
                lab = rownames(deg0vs1_noRibo),
                x = 'avg_log2FC',
                y = 'p_val_adj')
#EnhancedVolcano(deg0vs1_noRibo,
                lab = rownames(deg0vs1_noRibo),
                x = 'avg_log2FC',
                y = 'p_val_adj')
saveRDS(deg0vs1_noRibo, file = "deg0vs1_noRibo.rds")

##cluster1_2/0
deg0vs1_2 <- FindMarkers(microglia.combined_Umap, ident.1 = c(1,2), ident.2 =0, min.pct = 0.25)
head(deg0vs1_2)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
deg0vs1_2_noRibo<-deg0vs1_2[!grepl("^Rp[sl]",rownames(deg0vs1_2), ignore.case = F),]

deg0vs1_2_noRibo[order(deg0vs1_2_noRibo$p_val),]
library(EnhancedVolcano) 
EnhancedVolcano(deg0vs1_2_noRibo,ylim = c(0, 400),
                lab = rownames(deg0vs1_2_noRibo),
                x = 'avg_log2FC',
                y = 'p_val_adj')

saveRDS(deg0vs1_2_noRibo, file = "deg0vs1_2_noRibo.rds")

#cg_cluster124vs0.markers=cluster124vs0.markers[abs(cluster124vs0.markers$avg_log2FC) >0.4,] #log2FC?缍ㄖ?
#dim(cg_cluster124vs0.markers)
#cg_cluster124vs0.markers=cg_cluster124vs0.markers[order(cg_cluster124vs0.markers$avg_log2FC),]
#save(microglia1_2vs0.markers,file = 'microglia1_2vs0.markers.Rdata')


##
FeaturePlot(microglia.combined_Umap, features = c("S100a8" ,"S100a9"), min.cutoff = "q9",split.by = "seurat_clusters")
FeaturePlot(microglia.combined_Umap, features = c("Ppia" ,"Ptma"), min.cutoff = "q9",split.by = "seurat_clusters")
FeaturePlot(microglia.combined_Umap, features = c("Cst7" ,"Apoe"), min.cutoff = "q9",split.by = "seurat_clusters")
FeaturePlot(microglia.combined_Umap, features = c("Cst7" ,"Apoe"), min.cutoff = "q9",split.by = "orig.ident")
FeaturePlot(microglia.combined_Umap, features = c("Lpl" ,"Ccl2"), min.cutoff = "q9",split.by = "seurat_clusters")
top2<-microglia.combined.markers_noRibo%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
top2$gene
col<-paletteer_c("scico::berlin", 10)[c(1,3,4,9,5,2,6,8,10)]
DotPlot(microglia.combined_Umap,  features = c("Ppia","Fau","Ftl1","Eef1a1","Tpt1","Hmgb1","Pfn1", "Gapdh","Atp5g2", "Bin2"), cols = col, dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()

FeaturePlot(microglia.combined_Umap, features = c("Ifit2" ), min.cutoff = "q9",split.by = "orig.ident")



#下面代码没有跑
#杩欓噷鏄痬ergedata锛屾暟鎹泦闂村樊寮傚緢澶э紝鎵�浠ヤ笉鑳界洿鎺ョ畝鍗曡皟鐢‵indAllMarkers锛燂紵
#Macrophages:Cluster3/8
FeaturePlot(microglia.combined_Umap, features = c("Cd74","H2-Aa","H2-Ab1","H2-Eb1"), min.cutoff = "q9",split.by = "seurat_clusters")
#Macrophages:Cluster3,BAM亚群
FeaturePlot(microglia.combined_Umap, features = c("Pf4", "Dab2",  "F13a1" ), min.cutoff = "q9",split.by = "orig.ident")
FeaturePlot(microglia.combined_Umap, features = c("H2-Ab1", "Dab2",  "F13a1","Cd74","H2-Aa" ), min.cutoff = "q9",split.by = "seurat_clusters")
FeaturePlot(microglia.combined_Umap, features = c("Pf4", "Dab2",  "F13a1","Cybb","Lyz2"), min.cutoff = "q9",split.by = "seurat_clusters")
ggsave("BAM_genes_from_Cluster3.png",
       width = 7,             # 宽
       height = 7,            # 高
       units = "in",          # 单位
       dpi = 300    )          # 分辨率DPI)
#Macrophages:Cluster8,
Cluster8.markers <- FindConservedMarkers(microglia.combined, ident.1 = 8, grouping.var = "orig.ident", verbose = FALSE)
rownames(Cluster8.markers)
FeaturePlot(microglia.combined, features = c(  "Mrc1","Cd163","Cd74"), min.cutoff = "q9",split.by = "seurat_clusters")
FeaturePlot(microglia.combined, features = c(  "Mrc1","Cd163","Cd74"), min.cutoff = "q9",split.by = "orig.ident")
ggsave("PeriMacrophage_genes_from_Cluster8.png",
       width = 7,             # 宽
       height = 7,            # 高
       units = "in",          # 单位
       dpi = 300    )          # 分辨率DPI)

#Cluster1
Cluster1.markers <- FindConservedMarkers(microglia.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
rownames(Cluster1.markers)
FeaturePlot(microglia.combined, features = c("Ppia","Rpl21","Rps29","Hmgb1"), min.cutoff = "q9",split.by = "seurat_clusters")

ggsave("Microglia_genes_from_Cluster012.png",
       width = 7,             # 宽
       height = 7,            # 高
       units = "in",          # 单位
       dpi = 300    )          # 分辨率DPI)

#转为RNA后重新保存
saveRDS(microglia.combined, file = "microglia.combined_Umap.rds")

