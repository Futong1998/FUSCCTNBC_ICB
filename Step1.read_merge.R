#-------------------------------------------------------
####----------------0.environment prepare---------------
rm(list = ls())
options(stringsAsFactors = F)
gc()
dir.create("Step1.Filter")
setwd("./Step1.Filter")
#------------------------------------------------------
####----------------1.package prepare------------------
library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(DropletUtils)
library(scDblFinder)

Pre_sc_data_path <- "D:/nut/lab/TNBC_neo_scRNAseq/data/cellranger/Neo_pre"
Pre_sample_list <-  list.files(path=Pre_sc_data_path)
Pre_data.dir_list <- paste(Pre_sc_data_path,Pre_sample_list,"outs/filtered_feature_bc_matrix",sep = "/")

Late_sc_data_path <- "D:/nut/lab/TNBC_neo_scRNAseq/data/cellranger/late"
Late_sample_list <-  list.files(path=Late_sc_data_path)
Late_data.dir_path <- paste(Late_sc_data_path,Late_sample_list,"outs/filtered_feature_bc_matrix",sep = "/")

sc_data_path_list <- c(Pre_data.dir_list,Post_data.dir_list,Late_data.dir_path)
sample_list <- c(Pre_sample_list,Post_sample_list,Late_sample_list)
  


sc_data <- list()
for (i in 1:length(sample_list)) {
  Data <- Read10X(data.dir = sc_data_path_list[i])
  colnames(Data) <- paste(sample_list[i], colnames(Data), sep = '_')
  bcrank <- barcodeRanks(Data)
  # Only showing unique points for plotting speed.
  uniq <- !duplicated(bcrank$rank)
  tiff(filename = paste(sample_list[i],"bottomleft.tiff"))
  plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
  abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
  legend("bottomleft", legend=c("Inflection", "Knee"),
         col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  graphics.off()
  Data <- scDblFinder(Data)
  Data <- Data[,Data$scDblFinder.class=="singlet"]
  sc_data[[i]] <- Data
  print(i)
}

save(sc_data,file = "Step1.rmdl_36sample_scRNA_data.Rdata")

raw <- CreateSeuratObject(counts = cbind(sc_data[[1]]@assays@data@listData[["counts"]],
                                         sc_data[[2]]@assays@data@listData[["counts"]],
                                         sc_data[[3]]@assays@data@listData[["counts"]],
                                         sc_data[[4]]@assays@data@listData[["counts"]],
                                         sc_data[[5]]@assays@data@listData[["counts"]],
                                         sc_data[[6]]@assays@data@listData[["counts"]],
                                         sc_data[[7]]@assays@data@listData[["counts"]],
                                         sc_data[[8]]@assays@data@listData[["counts"]],
                                         sc_data[[9]]@assays@data@listData[["counts"]],
                                         sc_data[[10]]@assays@data@listData[["counts"]],
                                         sc_data[[11]]@assays@data@listData[["counts"]],
                                         sc_data[[12]]@assays@data@listData[["counts"]],
                                         sc_data[[13]]@assays@data@listData[["counts"]],
                                         sc_data[[14]]@assays@data@listData[["counts"]],
                                         sc_data[[15]]@assays@data@listData[["counts"]],
                                         sc_data[[16]]@assays@data@listData[["counts"]],
                                         sc_data[[17]]@assays@data@listData[["counts"]],
                                         sc_data[[18]]@assays@data@listData[["counts"]],
                                         sc_data[[19]]@assays@data@listData[["counts"]],
                                         sc_data[[20]]@assays@data@listData[["counts"]],
                                         sc_data[[21]]@assays@data@listData[["counts"]],
                                         sc_data[[22]]@assays@data@listData[["counts"]],
                                         sc_data[[23]]@assays@data@listData[["counts"]],
                                         sc_data[[24]]@assays@data@listData[["counts"]],
                                         sc_data[[25]]@assays@data@listData[["counts"]],
                                         sc_data[[26]]@assays@data@listData[["counts"]],
                                         sc_data[[27]]@assays@data@listData[["counts"]],
                                         sc_data[[28]]@assays@data@listData[["counts"]],
                                         sc_data[[29]]@assays@data@listData[["counts"]],
                                         sc_data[[30]]@assays@data@listData[["counts"]],
                                         sc_data[[31]]@assays@data@listData[["counts"]],
                                         sc_data[[32]]@assays@data@listData[["counts"]],
                                         project = "FUSCCTNBC_"))

length(table(raw@meta.data$orig.ident))
sum(table(raw@meta.data$orig.ident))
table(raw@meta.data$orig.ident)


raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-") 

tiff(filename = paste(sampleID,"_percent.tiff"),width = 1200,height = 400)
VlnPlot(raw,group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
tiff(filename = paste(sampleID,"_percent_cor.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()
plot1 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size = 0)
plot2 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0)
tiff(filename = paste(sampleID,"_percent_cor_prrm.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()

tiff(filename = paste(sampleID,"_percent_dotrm.tiff"),width = 1200,height = 400)
VlnPlot(raw,group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

tiff(filename = paste(sampleID,"_hist_nCount_RNA.tiff"),width = 600,height = 400)
plot1 <- hist(raw@meta.data[["nCount_RNA"]])
dev.off()

tiff(filename = paste(sampleID,"_hist_nFeature_RNA.tiff"),width = 600,height = 400)
plot2 <- hist(raw@meta.data[["nFeature_RNA"]])
dev.off()
tiff(filename = paste(sampleID,"_hist_percent.mt.tiff"),width = 600,height = 400)
plot3 <- hist(raw@meta.data[["percent.mt"]])
dev.off()

res <- raw@meta.data %>% dplyr::group_by(orig.ident) %>% dplyr::summarise(percent = sum(percent.mt>20)/n())

raw <- subset(raw, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20 & nCount_RNA> 400& nCount_RNA< 120000 )  #important step1: 要根据自己数据分布调整nFeature_RNA和percent.mt上限。
table(raw$orig.ident)
write.csv(table(raw$orig.ident),file = "post_filetre_cell_number.csv")

tiff(filename = paste(sampleID,"_percent_1.tiff"),width = 800,height = 400)
VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

plot1 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw,group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tiff(filename = paste(sampleID,"_percent_cor_2.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()


tiff(filename = paste(sampleID,"_percent_dotrm_2.tiff"),width = 2400,height = 400)
VlnPlot(raw,group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

tiff(filename = paste(sampleID,"_hist_nCount_RNA_2.tiff"),width = 600,height = 400)
plot1 <- hist(raw@meta.data[["nCount_RNA"]])
dev.off()

tiff(filename = paste(sampleID,"_hist_nFeature_RNA_2.tiff"),width = 600,height = 400)
plot2 <- hist(raw@meta.data[["nFeature_RNA"]])
dev.off()
tiff(filename = paste(sampleID,"_hist_percent.mt_2.tiff"),width = 600,height = 400)
plot3 <- hist(raw@meta.data[["percent.mt"]])
dev.off()


save(raw,file=paste(sampleID,"mt20_PrePCA.Rdata",sep = "_"))
gc()
########################################################################################
load("D:/nut/lab/TNBC_neo_scRNAseq/all_sample/FUSCCTNBC_35_PrePCA.Rdata")
npcs = 100
raw <- NormalizeData(raw, normalization.method = "LogNormalize", scale.factor = 10000)%>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = raw@var.genes, npcs = npcs, verbose = T)


VizDimLoadings(raw, dims = 1:2, reduction = "pca")

tiff(filename = paste(sampleID,npcs,"_PCA_pcs.tiff"),width = 800,height = 400)
DimPlot(raw, reduction = "pca")
dev.off()
tiff(filename = paste(sampleID,npcs,"_PCA_heatmap.tiff"),width = 1800,height = 3200)
DimHeatmap(raw, dims = 1:npcs, cells = 500, balanced = TRUE)
dev.off()

tiff(filename = paste(sampleID,npcs,"_ElbowPlot.tiff"),width = 1000,height = 500)
ElbowPlot(raw,ndims = npcs)
dev.off()


raw <- RunUMAP(raw, reduction = "pca", dims = 1:100)
raw <- FindNeighbors(raw, reduction = "pca", dims = 1:100)
raw <- FindClusters(raw, resolution = 0.4)

tiff(filename = paste(sampleID,npcs,"_umap_整合前.tiff"),width = 1800, height = 900)
p1 <- DimPlot(raw, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(raw, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

raw  <- raw %>% RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(raw, 'harmony')

tiff(filename = paste(sampleID,npcs,"harmony.tiff"),width =1800, height = 900)
p1 <- DimPlot(object = raw, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = raw, features = "harmony_1", group.by = "orig.ident",  pt.size = .1)
p1 + p2
dev.off()

tiff(filename = paste(sampleID,npcs,"_harmony_pcs.tiff"),width = 800,height = 400)
DimPlot(raw, reduction = "harmony")
dev.off()
tiff(filename = paste(sampleID,npcs,"_harmony_heatmap.tiff"),width = 1600,height = 6400)
DimHeatmap(raw, dims = 1:100, cells = 500, balanced = TRUE, reduction = "harmony")
dev.off()
tiff(filename = paste(sampleID,npcs,"harmony_ElbowPlot.tiff"),width = 1000,height = 500)
ElbowPlot(raw,ndims = npcs,reduction = 'harmony')
dev.off()


hpcs = 75
raw <- raw %>%
  RunUMAP(reduction = "harmony", dims = 1:hpcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:hpcs)

tiff(filename = paste(sampleID,hpcs,"_upcas.tiff"),width =900, height = 900)
DimHeatmap(raw, dims = 1:50, cells = 500, balanced = TRUE,reduction = "harmony")
dev.off()

raw<- CellCycleScoring(object = raw, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

save(raw,file =paste(sampleID,hpcs,"afterHarmony.Rdata",sep = "_"))








