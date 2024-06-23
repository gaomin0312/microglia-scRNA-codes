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

#Ê×ÏÈ£¬µ¼Èëdata
#¿ÉÒÔÌø¹ýÕâÀïµÄCHECK cells£¬Ö±½ÓloadÒÑ´¦ÀíµÄÊý¾Ý¼¯¡®sub_All6m_cluster0124_SCENIC_input.Rdata¡¯Êý¾Ý¼¯cluster0124monocle£¬À´Ô´ÈçÏÂ£º
microglia.combined_Umap
sce=microglia.combined_Umap
sce=subset(microglia.combined_Umap, downsample =100)
#sce NICåˆ†æž
sce  #çº?600ä¸ªç»†èƒžå·¦å³ï¼Œå¯èƒ½éœ€è¦è¿ç®—å¾ˆä¹…å¾ˆä¹?
#ä¸Šä¸€æ­¥éª¤ä¸­æå–å‡ºæ¥çš„æ„Ÿå…´è¶£äºšç¾?---2ç¾¤monocyteï¼Œæ¯ç¾¤æå?50ä¸ªåŸºå› è¿›è¡ŒåŽç»­è½¬å½•å› å­åˆ†æžï¼ˆä¸»è¦æ˜¯å› ä¸ºè¿ç®—é€Ÿåº¦å’Œæ—¶é—´çš„é™åˆ¶ï¼?
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
# ä¿è¯cisTarget_databases æ–‡ä»¶å¤¹ä¸‹é¢æœ‰ä¸‹è½½å¥? çš„æ–‡ä»?
#scenicOptions <- initializeScenic(org="mgi", dbDir=db, nCores=4)
scenicOptions<- initializeScenic(org="mgi", dbDir=db, nCores=1) #for runSENIC1-4³ÌÐò
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
# æœ€è€—è´¹æ—¶é—´çš„å°±æ˜¯è¿™ä¸ªæ­¥éª?

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
# è¿™ä¸ªæ­¥éª¤ä¹Ÿå¾ˆè€—æ—¶
scenicOptions<- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
# å› ä¸ºèŽ«åå…¶å¦™çš„é”™è¯¯ï¼Œéœ€è¦æŠŠ å¤šçº¿ç¨‹é‡æ–°è®¾ç½®æˆä¸? 1 ä¸ªçº¿ç¨?
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
#×÷ÕßÍÆ¼öµÄÔËËã½á¹û±£´æ£¬¡°Jimmy£ºÊµ¼ÊÉÏÓÃ²»ÉÏ¡±
#export2loom(scenicOptions, exprMat)  #ÔËÐÐºóÊäÈëloom
saveRDS(scenicOptions, file="microglia_int_scenicOptions.Rds") 

#½á¹û·ÖÎöºÍ¿ÉÊÓ»¯
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
motifEnrichment<-as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))# è¿è¡Œ 
#¿ÉÊÓ»ùÒòµÄmotifÐòÁÐÌØÕ÷
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Stat1"]
viewMotifs(tableSubset)
#¼ÓÉÏ»îÐÔµ¥Ôª£¨regulon£©µÄÏÞ¶¨ºó£¬motifÊýÁ¿¼õÉÙ
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat1" & highConfAnnot==TRUE]
viewMotifs(tableSubset)


#¿ÉÒÔÊ¹ÓÃµÄ½á¹û£º
#ÆäÖÐoutputÎÄ¼þ¼Ð
#±¾À´¾ÍÒÑ¾­×Ô¶¯»æÖÆÁË´óÁ¿µÄÍ¼±í¹©Ê¹ÓÃ£¬¶øÍ¼±í¶ÔÓ¦µÄÊý¾Ý¾Í´æ´¢ÔÚ loomFile ÀïÃæ

#½á¹û»æÖÆÍ¼±ê
##¿ÉÊÓ»¯
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