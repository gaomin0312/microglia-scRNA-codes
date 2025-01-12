rm(list = ls()) 
memory.limit(1000000)
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

setwd("D:/0scRNA_data_analysis/merge_test_noERCCscale/microglia")
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)
library(ComplexHeatmap)
options(width=200)

#首先，导入data
#可以跳过这里的CHECK cells，直接load已处理的数据集‘sub_All6m_cluster0124_SCENIC_input.Rdata’数据集cluster0124monocle，来源如下：
microglia.combined_Umap
sce=microglia.combined_Umap
sce=subset(microglia.combined_Umap, downsample =100)
#sce NIC鍒嗘瀽
sce  #绾?600涓粏鑳炲乏鍙筹紝鍙兘闇�瑕佽繍绠楀緢涔呭緢涔?
#涓婁竴姝ラ涓彁鍙栧嚭鏉ョ殑鎰熷叴瓒ｄ簹缇?---2缇onocyte锛屾瘡缇ゆ彁鍙?50涓熀鍥犺繘琛屽悗缁浆褰曞洜瀛愬垎鏋愶紙涓昏鏄洜涓鸿繍绠楅�熷害鍜屾椂闂寸殑闄愬埗锛?
#sce =subset(sce ,downsample=50)
#sce 
table(Idents(sce ))
phe=sce @meta.data   
mat=sce @assays$RNA@counts

mat[1:4,1:4]
exprMat =as.matrix(mat) 
dim(exprMat)
exprMat[1:4,1:4] 
head(phe)


cellInfo <-  phe[,c('seurat_clusters','nCount_RNA' ,'nFeature_RNA' )]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
cellInfo$CellType=Idents(sce )
table(cellInfo$CellType)

### Initialize settings
# https://github.com/aertslab/sce NIC
# https://pysce nic.readthedocs.io/en/latest/
library(SCENIC)
# https://resources.aertslab.org/cistarget/

db='D:/0scRNA_data_analysis/cisTarget_databases'
list.files(db)
# 淇濊瘉cisTarget_databases 鏂囦欢澶逛笅闈㈡湁涓嬭浇濂? 鐨勬枃浠?
#scenicOptions <- initializeScenic(org="mgi", dbDir=db, nCores=4)
scenicOptions<- initializeScenic(org="mgi", dbDir=db, nCores=1) #for runSENIC1-4程序
scenicOptions
saveRDS(cellInfo, file="cellInfo.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
length(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
#runCorrelation(exprMat_filtered, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
saveRDS(exprMat_filtered_log, file="exprMat_filtered_log.Rds")
# 鏈�鑰楄垂鏃堕棿鐨勫氨鏄繖涓楠?

gc()
setwd("D:/0scRNA_data_analysis/merge_test_noERCCscale/microglia")
weightMat <-runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(weightMat, file="weightMat.Rds")
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings

#showConnections(all = TRUE)
closeAllConnections()
scenicOptions<- runSCENIC_1_coexNetwork2modules(scenicOptions)
# 杩欎釜姝ラ涔熷緢鑰楁椂
scenicOptions<- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
# 鍥犱负鑾悕鍏跺鐨勯敊璇紝闇�瑕佹妸 澶氱嚎绋嬮噸鏂拌缃垚涓? 1 涓嚎绋?
#install.packages("doSNOW")
#install.packages("doParallel") 
#install.packages("doMPI")
#https://blog.csdn.net/m0_50412712/article/details/123614276
library(doSNOW)
library(doParallel)
library(doMPI)
scenicOptions<- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions<- runSCENIC_4_aucell_binarize(scenicOptions)
# Binary regulon activity: 29 TF regulons x 93 cells.
# (34 regulons including 'extended' versions)
# 29 regulons are active in more than 1% (0.93) cells.
tsneAUC(scenicOptions, aucType="AUC") # choose settings
#作者推荐的运算结果保存，“Jimmy：实际上用不上”
#export2loom(scenicOptions, exprMat)  #运行后输入loom
saveRDS(scenicOptions, file="microglia_int_scenicOptions.Rds") 

#结果分析和可视化
#rm(list = ls()) 
library(Seurat) 
library(SCENIC)
library(doParallel)
scenicOptions=readRDS(file="microglia_int_scenicOptions.Rds")

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org
# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
motifEnrichment<-as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))# 杩愯 
#可视基因的motif序列特征
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Stat1"]
viewMotifs(tableSubset)
#加上活性单元（regulon）的限定后，motif数量减少
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat1" & highConfAnnot==TRUE]
viewMotifs(tableSubset)


#可以使用的结果：
#其中output文件夹
#本来就已经自动绘制了大量的图表供使用，而图表对应的数据就存储在 loomFile 里面

#结果绘制图标
##可视化
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
scenicOptions=readRDS(file="microglia_int_scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name="Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)