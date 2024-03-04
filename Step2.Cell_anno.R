
##########################################
########    载入整合分群后数据 ###########
##########################################
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
sampleID <- "FUSCCTNBC"
setwd("/Step2.all_cell_anno")
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


