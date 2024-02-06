
##########################################
########    载入整合分群后数据 ###########
##########################################
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
sampleID <- "FUSCCTNBC_36"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step2.all_cell_anno")
##########################################
#############  大群差异分析 ##############
##########################################
markers            <- FindAllMarkers(raw,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers        <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

tiff(filename = paste(sampleID,"_marker_log2FC_heatmap.tiff"),width = 800,height = 1200)
DoHeatmap(raw, features = top_markers$gene) + NoLegend()
dev.off()

write.csv(markers,file = 'all_markers.csv')
write.csv(top_markers,file = 'top10_markers.csv')

##########################################
#############  特定marker表达情况 ########
##########################################
TNK_marker         <- c("CD3D","CD2","CD3E","CD4","CD8A","FOXP3")
B_marker           <- c("CD79A","MZB1","MS4A1","JCHAIN")
Myloid_marker      <- c("ITGAM","LYZ","CD68","TYROBP","HLA-DRA","CD14","FCGR3A","TPSAB1","CD86","CD163","CD33")
DC_marker          <- c("ITGAX","FCGR3A","CD209","IL3RA","HLA-DRA")
pDC_marker         <- c("LILRA4","CXCR3","IRF3")
mast_marker        <- c("CPA3","TPSAB1","TPSB2")
Epithelial_marker  <- c("KRT19","CD24","SCGB2A2")
Endo_marker        <- c("CLDN5","FLT1","RAMP2")
CAF_maker          <- c("COL1A1","DCN","C1R")
prolifer_marker    <- c("MKI67","STMN1")

p <- FeaturePlot(raw,features = TNK_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("TNK_marker_feature_plot",".tiff",sep = ""),width = 6,height = 9)

p <- FeaturePlot(raw,features = B_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("B_marker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = Myloid_marker,raster = T,ncol = 3)
ggsave(plot= p ,filename = paste("Myloid_marker_feature_plot",".tiff",sep = ""),width = 9,height = 12)

p <- FeaturePlot(raw,features = pDC_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("pDC_marker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = mast_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("mast_marker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = Epithelial_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("Epithelial_marker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = Endo_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("Endo_marker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = CAF_maker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("CAF_maker_feature_plot",".tiff",sep = ""),width = 6,height = 6)

p <- FeaturePlot(raw,features = prolifer_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("prolifer_marker_feature_plot",".tiff",sep = ""),width = 6,height = 3)


p <- FeaturePlot(raw,features = DC_marker,raster = T,ncol = 2)
ggsave(plot= p ,filename = paste("DC_marker_feature_plot",".tiff",sep = ""),width = 6,height = 3)


for (i in 1:length(TNK_marker)) {
  p <-  VlnPlot(raw,features = TNK_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("TNK_marker",TNK_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(B_marker)) {
  p <-  VlnPlot(raw,features = B_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("B_marker",B_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(prolifer_marker)) {
  p <-  VlnPlot(raw,features = prolifer_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("prolifer_marker",prolifer_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(Myloid_marker)) {
  p <-  VlnPlot(raw,features = Myloid_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("Myloid_marker",Myloid_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(Epithelial_marker)) {
  p <-  VlnPlot(raw,features = Epithelial_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("Epithelial_marker",Epithelial_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(Endo_marker)) {
  p <-  VlnPlot(raw,features = Endo_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("Endo_marker",Endo_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(CAF_maker)) {
  p <-  VlnPlot(raw,features = CAF_maker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("CAF_maker",CAF_maker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(pDC_marker)) {
  p <-  VlnPlot(raw,features = pDC_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("pDC_marker",pDC_marker[i],".tiff",sep = ""),width = 20,height = 6)
}

for (i in 1:length(mast_marker)) {
  p <-  VlnPlot(raw,features = mast_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("mast_marker",mast_marker[i],".tiff",sep = ""),width = 20,height = 6)
}
DimPlot(raw,label = T)
DimPlot(raw,cells.highlight =  WhichCells(raw,idents = c(39)),label = T)
DimPlot(raw,cells.highlight =  WhichCells(raw,idents = c(9,10,17,38,39)),label = T)
table(raw$seurat_clusters)
VlnPlot(raw,features = c("IFITM2","FCGR3B"),pt.size = 0)

for (i in 1:length(DC_marker)) {
  p <-  VlnPlot(raw,features = DC_marker[i],pt.size = 0.1)
  ggsave(plot= p ,filename = paste("DC_marker",DC_marker[i],".tiff",sep = ""),width = 20,height = 6)
}


tiff(filename = paste(sampleID,"_percent_serate.tiff"),width = 2400,height = 400)
VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

tiff(filename = paste(sampleID,"_percent_serate.tiff"),width = 2400,height = 400)
VlnPlot(raw, features = "IGHG3", pt.size = 0)
dev.off()

#9 10 17 33 38 39 49,51,54
C51_MARKER <- FindMarkers(raw,ident.1 = 51)
a <- c(0,1,4,6,8,13,15,20,22,41,42,43,45,53,55,56,57,
       5,14,16,24,27,29,31,38,
       11,12,18,32,33,49,
       2,3,9,10,17,23,28,34,35,36,37,40,44,46,47,50,51,52,54,
       7,19,25,26,30,39,
       21,48)
a[order(a)]
FeaturePlot(raw,features = c("MKI67","EPCAM","KRT19","CD3D"),raster = T,order = T)
DotPlot(raw,features = c("MKI67","EPCAM","KRT19","CD3D"))
FeaturePlot(subset(raw,KRT19==0),features = c("MKI67","EPCAM","CD3D","KRT19"),raster=T)

Feature

#######################
##髓系
library(readxl)
cellmarker_data <- read_excel("D:/nut/lab/TNBC_neo_scRNAseq/参考信息/cellmarker.xlsx")
gene_NM <- cellmarker_data[3:33,10:11]
colnames(gene_NM) <- c("GeneSymbol","class")

#######################
# 提取特定基因矩阵##NM文章
gene_NM <- subset(gene_NM,gene_NM$GeneSymbol%in%rownames(raw@assays$RNA@data))
expr_fanjia = raw@assays$RNA@data[gene_NM$GeneSymbol,]
# 制作注释文件
meta = raw@meta.data$seurat_clusters
names(meta)=rownames(raw@meta.data)
meta = as.data.frame(meta)
colnames(meta)='subcluster'
# 将同一cluster细胞的表达计算一个均值
new_expr_fanjia =do.call(rbind,
                         lapply(unique(meta$subcluster), function(t){
                           rowMeans(expr_fanjia[,meta$subcluster==t])
                         })
)
# 绘制分割热图
# 准备矩阵数据
matrix_fanjia = as.matrix(new_expr_fanjia)
scale_matrix_fanjia = scale(matrix_fanjia) # scale函数默认是按列，所以以列名为基因的形式进行scale
scale_matrix_fanjia = t(scale_matrix_fanjia)
colnames(scale_matrix_fanjia) = unique(meta$subcluster)

# 绘制分割热图
library(ComplexHeatmap)
pdf(file = paste(sampleID,"Myeliod_NM_gene_Anno.pdf",sep = ""),height = 12,width = 20)
Heatmap(scale_matrix_fanjia,cluster_rows = F,split = gene_NM$class,cluster_columns = T )
dev.off()



##########################################
#############  保存数据 ##################
##########################################    
write.csv(raw@meta.data[,'RNA_snn_res.1.4',drop=F],file = '4sample_clu_res1.4.csv')
save(raw,file = '5sample_afterFindMarkers.Rdata')

