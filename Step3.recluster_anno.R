rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

subcell <- subset(raw,merge_Cell_marjor=="TNK")
sampleID <- "T_NK"
dir.create("T_NK")
setwd("./T_NK")

subcell <- subset(raw,merge_Cell_marjor=="Myeloid")
sampleID <- "Myeloid"
dir.create("Myeloid")
setwd("./Myeloid")


subcell <- subset(raw,merge_Cell_marjor=="B")
sampleID <- "B"
dir.create("B")
setwd("./B")

subcell_1 <- NormalizeData(subcell, normalization.method = "LogNormalize", scale.factor = 10000)%>% ##????量???荼?准??,LogNormalize???惴???A = log( 1 + ( UMIA ?? UMITotal ) ?? 10000 )
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%##?????????弑?????(2000????,???????畏???,??PCA?????卟慰???NM????
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = subcell_1@var.genes, npcs = 100, verbose = FALSE)##PPCA??维??目???? ?_??


tiff(filename = paste(sampleID,"_ElbowPlot.tiff"),width = 1000,height = 500)
ElbowPlot(subcell_1,ndims = 100)
dev.off()


tiff(filename = paste(sampleID,"pcas_heatmap.tiff"),width =900, height = 3600)
DimHeatmap(subcell_1, dims = 1:100, cells = 500, balanced = TRUE)
dev.off()


##########去??????_endithelium??目???俨?去??
library(harmony)
subcell_3  <- subcell_1 %>% RunHarmony("orig.ident", plot_convergence = TRUE)

########################################################################
subcell_2 <- subcell_3 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)

subcell_2 <- FindClusters(subcell_2, resolution = 0.4)



markers            <- FindAllMarkers(subcell_2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers        <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.csv(markers,file = paste(sampleID,'_all_markers.csv'))
write.csv(top_markers,file = paste(sampleID,'_top10_markers.csv'))


library(readxl)
library(ggpubr)
library(purrr)
library(dplyr)

table(subcell_2$seurat_clusters)
result <- data.frame(
  celltype = names(table(subcell_2$seurat_clusters)))
GroupID <- names(table(subcell_2@meta.data$orig.ident_anno))
for (i in 1:length(GroupID)) {
  Data <- subcell_2@meta.data[subcell_2@meta.data$orig.ident_anno==GroupID[i],]
  freq <- table(Data$seurat_clusters)
  percent = as.data.frame(round(freq/sum(freq)*100, 2))
  colnames(percent) <- c("celltype",GroupID[i])
  result <- left_join(result,percent,by="celltype")
}
Data <- result
write.csv(Data,file = paste(sampleID,"propotion.csv"))

##########Early vs late
Data[is.na(Data)] <- 0
tmp_1 <- as.data.frame(t(Data[,c(pre_early,pre_late)]))
tmp_1$Group <- c(rep("Early",length(pre_early)),rep("Advanced",length(pre_late)))
colnames(tmp_1) <- c(Data$celltype,"Group")
tmp_data_1 <- reshape2::melt(tmp_1, "Group",variable_name = "Cell_cluster")
tmp_data_1$variable <- factor(tmp_data_1$variable,
                              levels = c("TNK","B","Myeloid","Mesenchymal","Endothelial","Normalepi","Cycling cell","Cancer"))

p <-ggboxplot(tmp_data_1, x = "variable", y = "value",
              color = "Group",palette =c("#A9747C","#9FAEC6"),add = "jitter")+
  ylab("Fraction")+xlab("Cell cluster")
p
ggsave(filename = "Stage_majoe.pdf")

res_1 <- compare_means(value~Group,data = tmp_data_1,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res_1, label = "p.signif",x='variable',y.position = 80)
p
ggsave(filename = "Stage_majoe.pdf")
pre_early_mean <- rowMeans(Data[res_1$variable,pre_early])
pre_late_mean<- rowMeans(Data[res_1$variable,pre_late])
res_1$pre_early_mean <- pre_early_mean
res_1$pre_late_mean <- pre_late_mean
res_1$pro_logFC <- log2(res_1$pre_early_mean/res_1$pre_late_mean)


library(ggrepel)
Dat <- res_1
Dat$p
Dat$threshold = factor(ifelse(Dat$p <0.05  & abs(Dat$pro_logFC) >= 0.5, 
                              ifelse(Dat$pro_logFC>= 0.5 ,'Higher in Early','Higher in Advacned'),'NoSignifi'),
                       levels=c('Higher in Early','Higher in Advacned','NoSignifi'))
Dat$all_men <- rowMeans(Data[res_1$variable,c(pre_early,pre_late)])

p <- ggplot(Dat,aes(x=pro_logFC,y=-log10(p),color=threshold,size=all_men))+
  geom_point()+
  scale_color_manual(values=c("#A5B5B3","#EC6399","#808080"))+#确定点的颜色
  geom_text_repel(
    data = Dat[Dat$p <0.05  & abs(Dat$pro_logFC) >= 0.5,],
    aes(label = variable),
    size = 3,color="black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p
ggsave(plot = p,filename = paste("FUSCCSCRNA_cellPor","EarlyvsAdv.pdf",sep = "_"),height = 5,width = 6)

write.csv(res_1,file = "FUSCCSCRNA_cellPor_propotion_WT_test_Stage.csv")
export::graph2pdf(p,file = paste(sampleID,"propotion_Stage.pdf"),width=9,height=6)

a <- res_1$variable[res_1$p.signif!="ns"]
a <- c("pDC","Mono.CXCL10")
for (i in 1:length(a)) {
  tmp_data <- tmp_data_1[tmp_data_1$variable %in% c(a[i]),]
  p <-ggboxplot(tmp_data, x = "Group", y = "value",color = "Group",palette =c("#9FAEC6","#A9747C"),add = "jitter")+
    ylab("Fraction")+xlab(a[i])+
    stat_compare_means()
  p
  export::graph2pdf(p,file = paste(sampleID,a[i],"_propotion_Stage.pdf"),width=3,height=6)
}
#####################
Data[is.na(Data)] <-0
tmp_4 <- as.data.frame(t(Data[,c(Chem_pCR,Chem_nonpCR,Immu_pCR,Immu_nonpCR)]))
tmp_4$Group <- c(rep("Chem_pCR",length(Chem_pCR)),rep("Chem_nonpCR",length(Chem_nonpCR)),rep("Immu_pCR",length(Immu_pCR)),rep("Immu_nonpCR",length(Immu_nonpCR)))
colnames(tmp_4) <- c(Data$celltype,"Group")
tmp_data_4 <- reshape2::melt(tmp_4, "Group",variable_name = "Cell_cluster")
names(table(tmp_data_4$variable))

p <-ggboxplot(tmp_data_4, x = "variable", y = "value",
              color = "Group",palette =c("#92AC9E","#D39B7B","#D08794","#667BAC"),add = "jitter")+
  ylab("Fraction")+xlab("Cell cluster")#+theme(axis.text.x = element_text(angle = 90,hjust = 1))
p
ggsave(plot = p,filename = "major_response.pdf")
export::graph2pdf(p,file = paste(sampleID,"propotion_REs_1.pdf"),width=16,height=6)


