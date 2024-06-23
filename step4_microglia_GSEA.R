rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

setwd("D:/0scRNA_data_analysis/merge_test_noERCCscale/sub_6m")


DefaultAssay(microglia.combined_Umap) <- "RNA"
library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
Loading required package: Rcpp

#msigdbr()函数可以得到想要的基因集，其中category参数可以指定基因集的类别

#小鼠差异基因，转大写
deg1_2_noRibo=microglia.combined.markers_noRibo[microglia.combined.markers_noRibo$cluster%in%c("1","2"),] #从FindAllMarker提取1/2而来
deg=deg1_2_noRibo
geneList= deg$avg_log2FC 
names(geneList)= toupper(rownames(deg))
names(geneList)
geneList=sort(geneList,decreasing = T)
geneList
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
#选择gmt文件（MigDB中的全部基因集）
gmtfile ='D:/0scRNA_data_analysis/merge_test_noERCCscale/sub_6m/MSigDB_gene_set/c7.all.v7.5.1.symbols.gmt'
# 31120 个基因集

#GSEA分析
library(GSEABase) # BiocManager::install('GSEABase')
geneset <- read.gmt( gmtfile )  
length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'microglia1_2_C7_gsea_results_df.csv')

library(enrichplot)
gseaplot2(egmt,geneSetID = 'GSE10325_BCELL_VS_MYELOID_DN', pvalue_table=TRUE)
gseaplot2(egmt,geneSetID = 'HALLMARK_MYC_TARGETS_V2',pvalue_table=TRUE) 
gseaplot2(egmt,geneSetID = 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',pvalue_table=T) 


#分析结果批量出图
down_kegg<-gsea_results_df[gsea_results_df$pvalue<0.05 & gsea_results_df$enrichmentScore < -0.3,];down_kegg$group=-1
up_kegg<-gsea_results_df[gsea_results_df$pvalue<0.05 & gsea_results_df$enrichmentScore > 0.3,];up_kegg$group=1


pro='test'
library(enrichplot) 
lapply(1:nrow(down_kegg), function(i){ 
  gseaplot2(egmt,down_kegg$ID[i],
            title=down_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_down_kegg_',
                gsub('/','-',down_kegg$Description[i])
                ,'.pdf'))
})
lapply(1:nrow(up_kegg), function(i){ 
  gseaplot2(egmt,up_kegg$ID[i],
            title=up_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_up_kegg_',
                gsub('/','-',up_kegg$Description[i]),
                '.pdf'))
})

