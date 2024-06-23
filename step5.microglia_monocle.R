rm(list = ls())
library(Seurat)
load(file = '...')

# 因为 monocle包是 使用 CellDataSet 对象，所以重新构�?
library(monocle)

# scater 
sample_ann <- sub.Flox.microglia_Umap@meta.data  
sample_ann$celltype=Idents(sub.Flox.microglia_Umap)
head(sample_ann)
#rownames(sample_ann)=sample_ann[,1]
gene_ann <- data.frame(
  gene_short_name = rownames(sub.Flox.microglia_Umap@assays$RNA) , 
  row.names =  rownames(sub.Flox.microglia_Umap@assays$RNA) 
)
head(gene_ann)

pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
ct=as.data.frame(sub.Flox.microglia_Umap@assays$RNA@counts)
ct[1:4,1:4]

sc_cds <- newCellDataSet(
  as.matrix(ct), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# 这个 CellDataSet 对象一定要认识清楚，务必花两个小时去摸索它�?


# 接下来仅仅是  monocle的标准流程而已
library(monocle)
sc_cds
sc_cds <- detectGenes(sc_cds, min_expr = 1) 
# 数值可以自行摸�?
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
sc_cds

cds <- sc_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) 
# plyr 

# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类�?
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1)
unsup_clustering_genes
cds <- setOrderingFilter(cds, 
                         unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, 
                           return_all = F) # norm_method='log'


# 其中 num_dim 参数选择基于上面的PCA�?
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 6) 
plot_cell_clusters(cds, 1, 2 )
table(pData(cds)$Cluster) 
colnames(pData(cds)) 
table(pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$seurat_clusters)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )

# 可以看到 monocle 给细胞重新定义了亚群，亚群数量是自己选择�?
# 整体来说，monocle和seurat 各自独立流程定义的亚群的一致性还不错�?

# 只是跑流程而已
save(cds,file = 'sub.Flox.microglia.monocle.input_cds.Rdata')

# 构建对象，seurat，monocle，scater
# monocle标准流程，降维聚类分�? 
rm(list = ls()) 
library(monocle)

load('microglia.monocle.input_cds.Rdat')
### 然后查看monocle ### 
cds 
# 接下来很重要，到底是看哪个性状的轨�?
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$celltype)
plot_cell_clusters(cds, 1, 2 )

## 我们这里并不能使�? monocle的分�?
# 还是依据前面�? seurat分群, 也就是说前面的代码仅仅是流程而已，我们没有使用那些结果哦

# 其实取决于自己真实的生物学意�?
pData(cds)$Cluster=pData(cds)$celltype
table(pData(cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
#  挑选差异最显著的基因可视化
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

# 前面是找差异基因，后面是做拟时序分析

# 第一�?: 挑选合适的基因. 有多个方法，例如提供已知的基因集�?
# 这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 
# 第二�?: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一�?, 其中“DDRTree”是Monocle2使用的默认方�?
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# 第三�?: 对细胞进行排�?
cds <- orderCells(cds)
# 最后两个可视化函数，对结果进行可视�?
plot_cell_trajectory(cds, color_by = "Cluster")  



length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster") 
#手动save png图片

phe=pData(cds)
boxplot(phe$Pseudotime,phe$Cluster)

# https://davetang.org/muse/2017/10/01/getting-started-monocle/
# 前面根据差异基因，推断好了拟时序，也就是说把差异基因动态化�?

# 后面就可以具体推断哪些基因随着拟时序如何的变化
my_cds_subset=cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1 )
# 不知道为什么在Mac电脑无法开启并行计算了 ，不过我测试了在Windows 电脑设置cores = 4是可以的
# 如果你是Mac电脑，自己修�? cores = 1 即可 
head(my_pseudotime_de)
save( my_cds_subset,my_pseudotime_de,
      file = 'microglia_monocle_output.Rdata')


#根据上述结果，更深进一步绘�?
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(gplots)
library(ggplot2)
library(monocle)
library(dplyr)
load(file = 'output_of_monocle.Rdata')
cds=my_cds_subset
phe=pData(cds)
colnames(phe)
library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm() 
p1

plot_cell_trajectory(cds, color_by = "celltype")  

p2=plot_cell_trajectory(cds, color_by = "Pseudotime")  
p2

p3=plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()
p3

library(patchwork)
p1+p2/p3

phe=pData(cds)
head(phe)
table(phe$State,phe$Cluster) 

library(dplyr)
my_pseudotime_de %>% arrange(qval) %>% head() 
# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene=my_pseudotime_gene[,1]
my_pseudotime_gene
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])+ scale_color_npg()
#已手动保存图片sub_All6m_Pseudotime_Cluster0124

plot_genes_jitter(my_cds_subset[my_pseudotime_gene,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )+ scale_color_nejm()
#已手动保存图片sub_All6m_Pseudotime_Cluster0124_3

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster
colnames(pData(my_cds_subset))
table(pData(my_cds_subset)$Cluster,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('celltype','State','Pseudotime')]
head(ac)
# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是�? Pseudotime 好了�?100列，50行的矩阵

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 # num_clusters = 2, 
                                                 # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster

print(my_pseudotime_cluster)
dev.off()


# 这个步骤超级耗费时间
# 不知道为什么在Mac电脑无法开启并行计算了 
# 不过我测试了在Windows 电脑设置cores = 4是可以的
# 如果你是Mac电脑，自己修�? cores = 1 即可 
BEAM_branch1 <- BEAM(my_cds_subset, branch_point = 1, cores = 4)
BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]
BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
head(BEAM_branch1) 

save(BEAM_branch1,file = 'sub_All6m_BEAM_res.Rdata')




# 使用全部的基因进行绘�? 
BEAM_res = BEAM_branch1
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)
#手动存图
pdf('monocle_BEAM_branch1_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()

#BEAM_res = BEAM_branch2
my_branched_heatmap <- plot_genes_branched_heatmap(
  my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
  branch_point = 1,
  num_clusters = 4, 
  use_gene_short_name = TRUE,
  show_rownames = F,
  return_heatmap = TRUE)

pdf('monocle_BEAM_branch2_heatmap.pdf')
print(my_branched_heatmap$ph)
dev.off()



head(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster) 
my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)

head(my_row[my_row$cluster == 3,'gene']) 

my_gene <- row.names(subset(fData(my_cds_subset),
                            gene_short_name %in% head(my_row[my_row$cluster == 1,'gene'])))
my_gene
# plot genes that are expressed in a branch dependent manner
plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 1,
                               ncol = 1)

#plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
branch_point = 2,
ncol = 1)
# 后面的批量绘图，意义不大 
names(pData(my_cds_subset))
head(pData(my_cds_subset))

plot_genes_jitter(my_cds_subset[gene_to_cluster,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow=  10,
                  ncol = NULL )

ggsave('monocle_top50_subCluster.pdf',height = 42)

plot_genes_in_pseudotime(my_cds_subset[head(gene_to_cluster,25),])
ggsave('monocle_top50_pseudotime.pdf',height = 49)


write.csv(my_pseudotime_de,file = 'microglia_pseudotime_de.csv')






