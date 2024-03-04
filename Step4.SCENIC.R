##############
#####SCENIC


library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(foreach)
rm(list=ls())
gc()



##准备细胞meta信息
dir.create("int")

raw <- subcell_2
cellInfo <- data.frame(raw@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="T_cell_recluster_anno")] <- "cluster"
cellInfo <- cellInfo[,c("sample","cluster")]


saveRDS(cellInfo, file="int/cellInfo_raw.Rds")
saveRDS(raw, "raw.Rds")

exprMat <- as.matrix(raw@assays$RNA@counts)

mydbDIR <- "/SCENIC/cis_download"
mydbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
           "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")


names(mydbs) <- c("500bp", "10kb")
library(parallel)
detectCores()
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "endo")
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
table(rawsub$B_cell_anno_V3)
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
save(exprMat_filtered_log,file = "exprMat_filtered_log.Rdata")

runGenie3(exprMat_filtered_log, scenicOptions,resumePreviousRun=T)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=20,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "endo")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
table(motifEnrichment_selfMotifs_wGenes$geneSet)
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="IRF1"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="STAT1" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "cluster"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
write.csv(rss,file = "rss.csv")


AUCmatrix  <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
rawauc <- AddMetaData(rawsub, AUCmatrix)
rawauc@assays$integrated <- NULL
saveRDS(rawauc,'rawauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
rawbin <- AddMetaData(rawsub, BINmatrix)
rawbin@assays$integrated <- NULL
saveRDS(rawbin, 'rawbin.rds')


######----------------------
##pheatmap
library(pheatmap)
library(dplyr)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(rawsub@meta.data,select = 'T_cell_anno')
celltype_or <-arrange(celltype,celltype$T_cell_anno)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
celltype_or_2 <- subset(celltype_or,rownames(celltype_or)%in%colnames(BINmatrix))
my.regulons <- c(colnames(rawauc@meta.data)[18:35])
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,rownames(celltype_or)]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,rownames(celltype_or_2)]

pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype_or,cluster_cols = F,use_raster=F,breaks = bk,
filename = 'scenic_seurat/myAUCmatrix_heatmap.pdf',
width = 6, height = 5)
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype_or_2,cluster_cols = G,use_raster=F,
filename = 'scenic_seurat/myBINmatrix_heatmap.pdf',
#color = colorRampPalette(colors = c("white","black"))(100),
width = 6, height = 5)
pheatmap(myBINmatrix,use_raster==T)
min(myBINmatrix)
max(myBINmatrix)





