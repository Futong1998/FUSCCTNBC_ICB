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


Post_sc_data_path <- "D:/nut/lab/TNBC_neo_scRNAseq/data/cellranger/Neo_post"
Post_sample_list <-  list.files(path=Post_sc_data_path)
Post_data.dir_list <- paste(Post_sc_data_path,Post_sample_list,"outs/filtered_feature_bc_matrix",sep = "/")

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
                                         sc_data[[33]]@assays@data@listData[["counts"]],
                                         sc_data[[34]]@assays@data@listData[["counts"]],
                                         sc_data[[35]]@assays@data@listData[["counts"]],
                                         sc_data[[36]]@assays@data@listData[["counts"]]),
                                         project = "FUSCCTNBC_36")

length(table(raw@meta.data$orig.ident))
sum(table(raw@meta.data$orig.ident))
table(raw@meta.data$orig.ident)

Late_sample <- c("CJL","S001","WFJ","XMG","ZFJ","ZQ")

anno_data <- data.frame(
  row.names = rownames(raw@meta.data),
  orig.ident_anno = as.character(raw$orig.ident),
  location =  rep("Breast",length(raw$orig.ident))
)

Post_sample_list
anno_data$orig.ident_anno[stringr::str_detect(rownames(anno_data),"CL_post")] <- "CL_post"
anno_data$orig.ident_anno[stringr::str_detect(rownames(anno_data),"CXR_post")] <- "CXR_post"
anno_data$orig.ident_anno[stringr::str_detect(rownames(anno_data),"SW_post_M")] <- "SW_post_M"
anno_data$orig.ident_anno[stringr::str_detect(rownames(anno_data),"SW_post_T")] <- "SW_post_T"
table(anno_data$orig.ident_anno)
length(table(anno_data$orig.ident_anno))

anno_data$Time <- ifelse(anno_data$orig.ident_anno%in%Pre_sample_list,"Pre_neo","Post_neo")
anno_data$Time [anno_data$orig.ident_anno%in%Late_sample] <- "Advanced"
table(anno_data$Time)

anno_data$location[anno_data$orig.ident_anno=="SW_post_M"|anno_data$orig.ident_anno=="WFJ"|anno_data$orig.ident_anno=="ZQ"] <- "ALN"
anno_data$location[anno_data$orig.ident_anno=="S001"|anno_data$orig.ident_anno=="ZFJ"] <- "Liver"
anno_data$location[anno_data$orig.ident_anno=="CJL"] <- "Chest"

table(anno_data$location)

anno_data$Stage <- ifelse(anno_data$orig.ident_anno%in%Late_sample,"Adcvanced","Early")
table(anno_data$Stage)


raw <- AddMetaData(raw,anno_data,col.name = colnames(anno_data))


save(raw,file = "35sample_Seurat_orig.Rdata")


table(anno_data$location)

##ȥ������������ͷ��ڴ�
a <- ls()
rm(list=a[which(a != 'raw')])
gc()
sampleID <- "FUSCCTNBC_36"
write.csv(table(raw$orig.ident_anno),file = "pre_filetre_cell_number.csv")
#load('27sample_Seurat_orig.Rdata')
###Seurat����
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-") #���gene�ļ��������������Сдmt���޸ĳ�Сдmt

##չʾ����������ٷֱ�
tiff(filename = paste(sampleID,"_percent.tiff"),width = 1200,height = 400)
VlnPlot(raw,group.by = "orig.ident_anno", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
tiff(filename = paste(sampleID,"_percent_cor.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()
plot1 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size = 0)
plot2 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0)
tiff(filename = paste(sampleID,"_percent_cor_prrm.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()

tiff(filename = paste(sampleID,"_percent_dotrm.tiff"),width = 1200,height = 400)
VlnPlot(raw,group.by = "orig.ident_anno", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
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

res <- raw@meta.data %>% dplyr::group_by(orig.ident_anno) %>% dplyr::summarise(percent = sum(percent.mt>20)/n())
# ���ݻ����������������������й���

raw <- subset(raw, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20 & nCount_RNA> 400& nCount_RNA< 120000 )  #important step1: Ҫ�����Լ����ݷֲ�����nFeature_RNA��percent.mt���ޡ�
table(raw$orig.ident_anno)
write.csv(table(raw$orig.ident_anno),file = "post_filetre_cell_number.csv")

tiff(filename = paste(sampleID,"_percent_1.tiff"),width = 800,height = 400)
VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

plot1 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw,group.by = "orig.ident_anno", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tiff(filename = paste(sampleID,"_percent_cor_2.tiff"),width = 800,height = 400)
CombinePlots(plots = list(plot1, plot2))
dev.off()


tiff(filename = paste(sampleID,"_percent_dotrm_2.tiff"),width = 2400,height = 400)
VlnPlot(raw,group.by = "orig.ident_anno", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
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

mean(raw$nCount_RNA)
mean(raw$nFeature_RNA)
mean(raw$percent.mt)

save(raw,file=paste(sampleID,"mt20_PrePCA.Rdata",sep = "_"))
gc()
########################################################################################
load("D:/nut/lab/TNBC_neo_scRNAseq/all_sample/FUSCCTNBC_35_PrePCA.Rdata")
npcs = 100
raw <- NormalizeData(raw, normalization.method = "LogNormalize", scale.factor = 10000)%>% ##���������ݱ�׼��,LogNormalize���㷨��A = log( 1 + ( UMIA �� UMITotal ) �� 10000 )
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%##��������߱����(2000����,�������η���,��PCA����߲ο���NM����
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = raw@var.genes, npcs = npcs, verbose = T)##PPCA��ά��Ŀ���� �_��

#####ѡ��ά����
##չʾ���ɷֻ����ֵ
VizDimLoadings(raw, dims = 1:2, reduction = "pca")
##����pcaɢ��ͼ
tiff(filename = paste(sampleID,npcs,"_PCA_pcs.tiff"),width = 800,height = 400)
DimPlot(raw, reduction = "pca")
dev.off()
tiff(filename = paste(sampleID,npcs,"_PCA_heatmap.tiff"),width = 1800,height = 3200)
DimHeatmap(raw, dims = 1:npcs, cells = 500, balanced = TRUE)
dev.off()

#####################################ֻ�ܳʬFǰ20��
# raw <- JackStraw(raw, num.replicate = 100)
# raw <- ScoreJackStraw(raw, dims = 1:20)
# tiff(filename = paste(sampleID,"_JackStrawPlot.tiff"),width = 1000,height = 500)
# JackStrawPlot(raw, dims = 1:20)
# dev.off()
#####################

tiff(filename = paste(sampleID,npcs,"_ElbowPlot.tiff"),width = 1000,height = 500)
ElbowPlot(raw,ndims = npcs)
dev.off()


###PCA��Ŀ����ǰ���ͼ����ȷ��������ͼѡ��ǰ50
raw <- RunUMAP(raw, reduction = "pca", dims = 1:100)
raw <- FindNeighbors(raw, reduction = "pca", dims = 1:100)
raw <- FindClusters(raw, resolution = 0.4)

tiff(filename = paste(sampleID,npcs,"_umap_����ǰ.tiff"),width = 1800, height = 900)
p1 <- DimPlot(raw, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(raw, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

save(raw,file=paste(sampleID,npcs,"Preharmony.Rdata",sep = "_"))
## ����
raw  <- raw %>% RunHarmony("orig.ident_anno", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(raw, 'harmony')

##����֮֮��� �ص���
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


##���γ���
hpcs = 75
raw <- raw %>%
  RunUMAP(reduction = "harmony", dims = 1:hpcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:hpcs)

tiff(filename = paste(sampleID,hpcs,"_upcas_���Ϻ�.tiff"),width =900, height = 900)
DimHeatmap(raw, dims = 1:50, cells = 500, balanced = TRUE,reduction = "harmony")
dev.off()


raw <- FindClusters(raw, resolution = 0.2)
raw <- FindClusters(raw, resolution = 0.4)
raw <- FindClusters(raw, resolution = 0.6)
raw <- FindClusters(raw, resolution = 0.8)
raw <- FindClusters(raw, resolution = 1.0)
raw <- FindClusters(raw, resolution = 1.2)
raw <- FindClusters(raw, resolution = 1.4)


tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_patient.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "orig.ident_anno")
dev.off()

# tiff(filename = paste(sampleID,"_umap_���Ϻ�_Group.tiff"),width =1800, height = 1900)
# DimPlot(raw, reduction = "umap", group.by = "Group")
# dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_0.2.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.0.2", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_0.4.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.0.4", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_0.6.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.0.6", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_0.8.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_1.0.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.1.0", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_1.2.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.1.2", label = TRUE, repel = TRUE)
dev.off()

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_clus_1.4.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "RNA_snn_res.1.4", label = TRUE, repel = TRUE)
dev.off()


#Idents(raw) <- 'RNA_snn_res.0.4'


##--ϸ������
raw<- CellCycleScoring(object = raw, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

tiff(filename = paste(sampleID,hpcs,"_umap_���Ϻ�_cellcycle.tiff"),width =900, height = 900)
DimPlot(raw, reduction = "umap", group.by = "Phase", label = TRUE, repel = TRUE)
dev.off()

# saveRDS(raw, file = "6sample.rds")
# save(raw,file="6sample.Robj")##����Ҫ����һ���켣������Ҫ�õ�
save(raw,file =paste(sampleID,hpcs,"afterHarmony.Rdata",sep = "_"))







