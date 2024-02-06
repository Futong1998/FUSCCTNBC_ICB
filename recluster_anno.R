rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

sampleID <- "T_NK"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res0.4")
load("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/T_NK_afterharmony.Rdata")

sampleID <- "T_NK_RES0.7"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res.7")

sampleID <- "Myeloid"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/Myeloid")
########################################################################
subcell_2 <- subcell_3 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20)

subcell_2 <- FindClusters(subcell_2, resolution = 0.4)
subcell_2 <- FindClusters(subcell_2, resolution = 0.5)
subcell_2 <- FindClusters(subcell_2, resolution = 0.7)
subcell_2 <- FindClusters(subcell_2, resolution = 0.8)
subcell_2 <- FindClusters(subcell_2, resolution = 1)
subcell_2 <- FindClusters(subcell_2, resolution = 1.2)
subcell_2 <- FindClusters(subcell_2, resolution = 1.4)

sampleID <- "T_NK_res0.4"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res0.4")
Idents(subcell_2) <- subcell_2$RNA_snn_res.0.4
subcell_2$seurat_clusters <- subcell_2$RNA_snn_res.0.4

sampleID <- "T_NK_RES0.7"
setwd("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res.7")
Idents(subcell_2) <- subcell_2$RNA_snn_res.0.7
subcell_2$seurat_clusters <- subcell_2$RNA_snn_res.0.7

markers            <- FindAllMarkers(subcell_2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers        <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.csv(markers,file = paste(sampleID,'_all_markers.csv'))
write.csv(top_markers,file = paste(sampleID,'_top10_markers.csv'))

save(subcell_2,file = paste(sampleID,"recluster.Rdata",sep = "_"))
p <- DimPlot(subcell_2)
ggsave(plot = p,filename = paste(sampleID,"_dimplot.pdf"),width = 7,height = 6)

p <- FeaturePlot(subcell_2,c('CD4',"CD8A","GZMB","CD3E"))
ggsave(plot = p,filename = paste(sampleID,"_T_submarker_feat.pdf"),width = 6,height = 6)

p <- VlnPlot(subcell_2,c('CD4',"CD8A"))
ggsave(plot = p,filename = paste(sampleID,"_T_submarker.pdf"),width = 12,height = 6)


DoHeatmap(subcell_2,features = top_markers$gene)
ggsave(filename = "Doheatmap_top5marker.pdf",width = 12,height = 8)
############################################################################################################
##细胞比例类型
library(readxl)
library(ggpubr)
library(purrr)
library(dplyr)
##Load_data
Patient_data <- read_excel("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Group_for_sc_Analysis.xlsx")
pre_early <- Patient_data$Sample[Patient_data$Time=="Pre_neo"]
post_early <- Patient_data$Sample[Patient_data$Time=="Post_neo"]
pre_late <- Patient_data$Sample[Patient_data$Time=="Advanced"]
pCR <- Patient_data$Sample[Patient_data$Response=="pCR"]
non_PCR <- Patient_data$Sample[Patient_data$Response=="nonpCR"]
Chem_pCR <-  Patient_data$Sample[Patient_data$Response=="pCR"&Patient_data$Group=="Chem"]
Chem_nonpCR<- Patient_data$Sample[Patient_data$Response=="nonpCR"&Patient_data$Group=="Chem"]
Immu_pCR<- Patient_data$Sample[Patient_data$Response=="pCR"&Patient_data$Group=="Immune"]
Immu_nonpCR<- Patient_data$Sample[Patient_data$Response=="nonpCR"&Patient_data$Group=="Immune"]


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
tmp_1 <- as.data.frame(t(Data[,c(pre_early,pre_late)]))
tmp_1$Group <- c(rep("Early",length(pre_early)),rep("Advanced",length(pre_late)))
colnames(tmp_1) <- c(Data$celltype,"Group")
tmp_data_1 <- reshape2::melt(tmp_1, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_1, x = "variable", y = "value",
              color = "Group",palette =c("#EC6399","#9E5F33"))+
  ylab("Fraction")+xlab("Cell cluster")
p

res_1 <- compare_means(value~Group,data = tmp_data_1,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res_1, label = "p.signif",x='variable',y.position = 40)
p

export::graph2pdf(p,file = paste(sampleID,"propotion_Stage.pdf"),width=12,height=6)

tmp_data <- tmp_data_1[tmp_data_1$variable %in% c(3,5,7,12),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#EC6399","#9E5F33"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 30)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_Stage_sig.pdf"),width=8,height=6)

#########pCR vs nonPCR(all) 
tmp_2 <- as.data.frame(t(Data[,c(pCR,non_PCR)]))
tmp_2$Group <- c(rep("pCR",length(pCR)),rep("non_PCR",length(non_PCR)))
colnames(tmp_2) <- c(Data$celltype,"Group")

tmp_data_2 <- reshape2::melt(tmp_2, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_2, x = "variable", y = "value",
              color = "Group",palette =c("#446EB3","#E65634"))+
  ylab("Fraction")+xlab("Cell cluster")
res_2 <- compare_means(value~Group,data = tmp_data_2,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res_2, label = "p.signif",x='variable',y.position = 40)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_pCR_all.pdf"),width=12,height=6)

tmp_data <- tmp_data_2[tmp_data_2$variable %in% c(8),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#446EB3","#E65634"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_pCR_sig.pdf"),width=4,height=6)


#########pCR vs nonPCR(Chem) 
tmp_3 <- as.data.frame(t(Data[,c(Chem_pCR,Chem_nonpCR)]))
tmp_3$Group <- c(rep("Chem_pCR",length(Chem_pCR)),rep("Chem_nonpCR",length(Chem_nonpCR)))
colnames(tmp_3) <- c(Data$celltype,"Group")
tmp_data_3 <- reshape2::melt(tmp_3, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_3, x = "variable", y = "value",
              color = "Group",palette =c("#96B8DB","#F58268"))+
  ylab("Fraction")+xlab("Cell cluster")
res_3 <- compare_means(value~Group,data = tmp_data_3,method = "wilcox.test",group.by = "variable")    
p <- p+stat_pvalue_manual(res_3, label = "p.signif",x='variable',y.position = 30)
p

export::graph2pdf(p,file = paste(sampleID,"propotion_Chem_pCR_all.pdf"),width=12,height=6)

#########pCR vs nonPCR(Immune) 
tmp_4 <- as.data.frame(t(Data[,c(Immu_pCR,Immu_nonpCR)]))
tmp_4$Group <- c(rep("Immu_pCR",length(Immu_pCR)),rep("Immu_nonpCR",length(Immu_nonpCR)))
colnames(tmp_4) <- c(Data$celltype,"Group")
tmp_data_4 <- reshape2::melt(tmp_4, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_4, x = "variable", y = "value",
              color = "Group",palette =c("#00A4E8","#F4979C") )+
  ylab("Fraction")+xlab("Cell cluster")
res_4 <- compare_means(value~Group,data = tmp_data_4,method = "wilcox.test",group.by = "variable")   
p <- p+stat_pvalue_manual(res_4, label = "p.signif",x='variable',y.position = 30)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_Immu_pCR.pdf"),width=12,height=6)

tmp_data <- tmp_data_4[tmp_data_4$variable %in% c(5),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#00A4E8","#F4979C"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_Immu_pCR_sig.pdf"),width=4,height=6)


#########pCR vs nonPCR(Immune) 
Data <- read.csv("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res.7/细胞比例_比较/T_NK_RES0.7 propotion.csv",row.names = 1)
Data <- read_excel("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/Myeloid/Myeloid/output_pc25res0.8/Myeloid propotion.xlsx")
Data <- read.csv("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step6.tumor_identi/cANCER propotion.csv",row.names = 1)

tmp_5 <- as.data.frame(t(Data[,c(Immu_pCR,Immu_nonpCR,Chem_pCR,Chem_nonpCR)]))
tmp_5$Group <- c(rep("Immu_pCR",length(Immu_pCR)),rep("Immu_nonpCR",length(Immu_nonpCR)),rep("Chem_pCR",length(Chem_pCR)),rep("Chem_nonpCR",length(Chem_nonpCR)))
colnames(tmp_5) <- c(Data$celltype,"Group")
tmp_data_5 <- reshape2::melt(tmp_5, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_5, x = "variable", y = "value",
              color = "Group",palette =c("#6365E1","#FF7B02","#00A4E8","#F4979C") )+
  ylab("Fraction")+xlab("Cell cluster")
res_5 <- compare_means(value~Group,data = tmp_data_5,method = "wilcox.test",group.by = "variable")   
p
export::graph2pdf(p,file = paste("Cancer_propotion_Immu_pCR.pdf",sep = ""),width=12,height=4)

export::graph2pdf(p,file = paste("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step3.recluster_V2/T_NK/res.7/细胞比例_比较/","propotion_Immu_pCR.pdf",sep = ""),width=8,height=4)

tmp_data <- tmp_data_5[tmp_data_5$variable %in% c(0,1,7),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#6365E1","#FF7B02","#00A4E8","#F4979C"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"propotion_Immu_pCR_sig.pdf"),width=6,height=6)
export::graph2ppt(p,file = paste(sampleID,"propotion_Immu_pCR_sig.ppt"),width=4,height=6)

# #######pre vs pOST NEO    
# tmp_5 <- as.data.frame(t(Data[,c(pre_early,post_early)]))
# tmp_5$Group <- c(rep("pre_early",length(pre_early)),rep("post_early",length(post_early)))
# colnames(tmp_5) <- c(Data$celltype,"Group")
# tmp_data_5 <- reshape2::melt(tmp_5, "Group",variable_name = "Cell_cluster")
# p <-ggboxplot(tmp_data_5, x = "variable", y = "value",
#               color = "Group",palette =c("#FEE290","#8B3B99"))+
#   ylab("Fraction")+xlab("Cell cluster")
# res_5 <- compare_means(value~Group,data = tmp_data_5,method = "wilcox.test",group.by = "variable")  
# p <- p+stat_pvalue_manual(res_5, label = "p.signif",x='variable',y.position = 60)
# p
# export::graph2pdf(p,file = paste(sampleID,"propotion_time.pdf"),width=12,height=6)
# 
# tmp_data <- tmp_data_2[tmp_data_2$variable %in% c(3),]
# p <-ggboxplot(tmp_data, x = "variable", y = "value",
#               color = "Group",palette =c("#FEE290","#8B3B99"))+
#   ylab("Fraction")+xlab("Cell cluster")
# res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
# p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 50)
# p
# export::graph2pdf(p,file = paste(sampleID,"propotion_time_select.pdf"),width=4,height=6)
# 

result_test <- list(res_1,res_2,res_3,res_4,res_5) %>% reduce(full_join, by = "variable")
write.csv(result_test,file = paste(sampleID,'subtype_compare.csv'))


#######--------绝对比例--------------------------------   
cell_pro <- major_pro[major_pro$celltype=="TNK",2:37]
Data_abs <- result[,2:37]
Data_abs_1 <- Data_abs
for (a in 1:ncol(Data_abs)) {
  Data_abs_1[,a] <- Data_abs[,a]*cell_pro[,a]/100
}
##########Early vs late
tmp_1 <- as.data.frame(t(Data_abs_1[,c(pre_early,pre_late)]))
tmp_1$Group <- c(rep("pre_early",length(pre_early)),rep("pre_late",length(pre_late)))
colnames(tmp_1) <- c(Data$celltype,"Group")
tmp_data_1 <- reshape2::melt(tmp_1, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_1, x = "variable", y = "value",
              color = "Group",palette =c("#EC6399","#9E5F33"))+
  ylab("Fraction")+xlab("Cell cluster")
res_1 <- compare_means(value~Group,data = tmp_data_1,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res_1, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Stage.pdf"),width=12,height=6)

tmp_data <- tmp_data_1[tmp_data_1$variable %in% c(3,4,5,7,12),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#EC6399","#9E5F33"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 15)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Stage_sig.pdf"),width=8,height=6)

#########pCR vs nonPCR(all) 
tmp_2 <- as.data.frame(t(Data_abs_1[,c(pCR,non_PCR)]))
tmp_2$Group <- c(rep("pCR",length(pCR)),rep("non_PCR",length(non_PCR)))
colnames(tmp_2) <- c(Data$celltype,"Group")

tmp_data_2 <- reshape2::melt(tmp_2, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_2, x = "variable", y = "value",
              color = "Group",palette =c("#446EB3","#E65634"))+
  ylab("Fraction")+xlab("Cell cluster")
res_2 <- compare_means(value~Group,data = tmp_data_2,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res_2, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_pCR.pdf"),width=12,height=6)

tmp_data <- tmp_data_2[tmp_data_2$variable %in% c(2),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#446EB3","#E65634")+
                ylab("Fraction")+xlab("Cell cluster"))
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 5)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_pCR_sig.pdf"),width=4,height=6)

#########pCR vs nonPCR(Chem) 
tmp_3 <- as.data.frame(t(Data_abs_1[,c(Chem_pCR,Chem_nonpCR)]))
tmp_3$Group <- c(rep("Chem_pCR",length(Chem_pCR)),rep("Chem_nonpCR",length(Chem_nonpCR)))
colnames(tmp_3) <- c(Data$celltype,"Group")
tmp_data_3 <- reshape2::melt(tmp_3, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_3, x = "variable", y = "value",
              color = "Group",palette =c("#96B8DB","#F58268"))+
  ylab("Fraction")+xlab("Cell cluster")
res_3 <- compare_means(value~Group,data = tmp_data_3,method = "wilcox.test",group.by = "variable")    
p <- p+stat_pvalue_manual(res_3, label = "p.signif",x='variable',y.position = 30)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Chem_pCR.pdf"),width=12,height=6)

#########pCR vs nonPCR(Immune) 
tmp_4 <- as.data.frame(t(Data_abs_1[,c(Immu_pCR,Immu_nonpCR)]))
tmp_4$Group <- c(rep("Immu_pCR",length(Immu_pCR)),rep("Immu_nonpCR",length(Immu_nonpCR)))
colnames(tmp_4) <- c(Data$celltype,"Group")
tmp_data_4 <- reshape2::melt(tmp_4, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_4, x = "variable", y = "value",
              color = "Group",palette =c("#00A4E8","#F4979C") )+
  ylab("Fraction")+xlab("Cell cluster")
res_4 <- compare_means(value~Group,data = tmp_data_4,method = "wilcox.test",group.by = "variable")   
p <- p+stat_pvalue_manual(res_4, label = "p.signif",x='variable',y.position = 20)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Immu_pCR.pdf"),width=12,height=6)

tmp_data <- tmp_data_4[tmp_data_4$variable %in% c(5),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#00A4E8","#F4979C"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 8)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Immu_pCR_sig.pdf"),width=4,height=6)



#######pre vs pOST NEO    
tmp_5 <- as.data.frame(t(Data_abs_1[,c(pre_early,post_early)]))
tmp_5$Group <- c(rep("pre_early",length(pre_early)),rep("post_early",length(post_early)))
colnames(tmp_5) <- c(Data$celltype,"Group")
tmp_data_5 <- reshape2::melt(tmp_5, "Group",variable_name = "Cell_cluster")
p <-ggboxplot(tmp_data_5, x = "variable", y = "value",
              color = "Group",palette =c("#FEE290","#8B3B99"))+
  ylab("Fraction")+xlab("Cell cluster")
res_5 <- compare_means(value~Group,data = tmp_data_5,method = "wilcox.test",group.by = "variable")  
p <- p+stat_pvalue_manual(res_5, label = "p.signif",x='variable',y.position = 30)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Time.pdf"),width=12,height=6)

tmp_data <- tmp_data_5[tmp_data_5$variable %in% c(5,10),]
p <-ggboxplot(tmp_data, x = "variable", y = "value",
              color = "Group",palette =c("#FEE290","#8B3B99"))+
  ylab("Fraction")+xlab("Cell cluster")
res <- compare_means(value~Group,data = tmp_data,method = "wilcox.test",group.by = "variable")
p <- p+stat_pvalue_manual(res, label = "p.signif",x='variable',y.position = 5)
p
export::graph2pdf(p,file = paste(sampleID,"abs_propotion_Time_sig.pdf"),width=4,height=6)

result_test_abs <- list(res_1,res_2,res_3,res_4,res_5) %>% reduce(full_join, by = "variable")
write.csv(result_test_abs,file = paste(sampleID,'abs_propotion_Time_sig.pdf'))

#########
##新辅助配对样本
#--------相对
dat <- as.data.frame(t(Data[,c("CL","CL_post","SW","SW_post_T","CXR","CXR_post")]))
dat$Group <- c("pre","post","pre","post","pre","post")
colnames(dat) <- c(Data$celltype,"Group")

order <- c("CL","CL","SW","SW","CXR","CXR")

dat_1 <- reshape2::melt(dat, "Group",variable_name = "Cell_cluster")
dat_1$patient <- rep(order,nrow(Data))
colnames(dat_1) <- c("Group","Cell","Propotion","Patient")
dat_2 <- dplyr::arrange(dat_1,dat_1$Group,dat_1$Patient)
dat_2$Group <- factor(dat_2$Group,levels = c("pre","post"))
p <- ggpaired(dat_2, x = "Group", y = "Propotion",
              color = "Group", palette = c("#8B3B99","#FEE290"),
              line.color = "gray", line.size = 0.4,
              facet.by = "Cell", short.panel.labs = F)+ 
  stat_compare_means(label = "p.format", paired = TRUE,label.y = 30)
p
ggsave(plot = p,filename = paste(sampleID,"Pre_after_neo_cellprop.pdf",sep = "_"),width = 12,height = 12)

#--------绝对
dat <- as.data.frame(t(Data_abs_1[,c("CL","CL_post","SW","SW_post_T","CXR","CXR_post")]))
dat$Group <- c("pre","post","pre","post","pre","post")
colnames(dat) <- c(Data$celltype,"Group")

order <- c("CL","CL","SW","SW","CXR","CXR")

dat_1 <- reshape2::melt(dat, "Group",variable_name = "Cell_cluster")
dat_1$patient <- rep(order,nrow(Data))
colnames(dat_1) <- c("Group","Cell","Propotion","Patient")
dat_2 <- dplyr::arrange(dat_1,dat_1$Group,dat_1$Patient)
dat_2$Group <- factor(dat_2$Group,levels = c("pre","post"))
p <- ggpaired(dat_2, x = "Group", y = "Propotion",
              color = "Group", palette = c("#8B3B99","#FEE290"),
              line.color = "gray", line.size = 0.4,
              facet.by = "Cell", short.panel.labs = F)+ 
  stat_compare_means(label = "p.format", paired = TRUE,label.y = 6)
p
ggsave(plot = p,filename = paste(sampleID,"Pre_after_neo_abs_cellprop.pdf",sep = "_"),width = 12,height = 12)



############################################################################################################
##细胞类型注释
expr_data = subcell_2@assays$RNA@data[top_markers$gene,]
# ????????????
meta = subcell_2$seurat_clusters
names(meta)=rownames(subcell_2@meta.data)
meta = as.data.frame(meta)
colnames(meta)='cell_type'
meta <- subset(meta,!is.na(meta$cell_type))
# ??????cluster??????????????????????
new_expr_data =do.call(rbind,
                       lapply(unique(meta$cell_type), function(t){
                         rowMeans(expr_data[,meta$cell_type==t])
                       })
)
# ????????

library(pheatmap)
matrix = as.matrix(new_expr_data)
scale_matrix = scale(matrix) # scale??????????????????????????????????????????scale
scale_matrix = t(scale_matrix)
colnames(scale_matrix) = as.character(unique(meta$cell_type))
scale_matrix <- as.data.frame(scale_matrix[,as.character(unique(top_markers$cluster))])
bk = unique(c(seq(-2,2, length=100)))
mycolors <- c(colorRampPalette(c("#74ADD1", "white"))(50), colorRampPalette(c("white", "#F46D43"))(50))

p <- pheatmap(as.matrix(scale_matrix),cluster_rows = F,cluster_cols = F, color = mycolors,border_color =NA,breaks = bk )+coor
p
export::graph2ppt(p,file=paste(sampleID,"_top5marker.pptx"),width = 6,height = 8)

RColorBrewer::brewer.pal()
p <- DotPlot(subcell_2,features = unique(top_markers$gene))+
  theme(axis.text.x = element_text(angle = 90))
p
ggsave(plot = p,filename = paste(sampleID,"Top5gene_marker.pdf"),height = 5,width = 12)

ggsave(plot = p,filename = paste(sampleID,"Top5gene_clu.pdf"),height = 8,width = 6)

p <- DoHeatmap(subcell_2,features = top_markers$gene,label = F,slot = "scale.data")
p
ggsave(plot = p,filename = paste(sampleID,"Top5gene_heatmap.pdf"),height = 8,width = 8)


# ----------- ???? ------------- ##
# ????????????
###Fanjia
gene_cell_type = c('CD3G','CD3D','CD3E','CD4','CD8A','CD8B','KLRF1','KLRD1')
gene_naive = c('TCF7','SELL','LEF1','CCR7')
gene_inhibitory = c('LAG3','TIGIT','PDCD1','HAVCR2','CTLA4')
gene_cytotoxic = c('IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7','IL17A','IL17F','LAMTOR3')
gene_co_stimulatory = c('CD28','TNFRSF14','ICOS','TNFRSF9')
gene_resident = c('CD69','RUNX3','NR4A1')
gene_TF = c('TBX21','ZNF683','ZEB2','ID2','EOMES','HOPX','HIF1A','TOX')
gene_Treg = c('FOXP3','IL2RA','IKZF2')
gene_fanjia = c(gene_cell_type,gene_naive,gene_inhibitory,gene_cytotoxic,gene_co_stimulatory,gene_resident,gene_TF,gene_Treg)

##SA
gene_SA_CD4_memory = c('ANXA1')
gene_SA_CD4_effector_memory = c('ANXA1','GZMA')
gene_SA_CD4_exhaustd = c('CXCL13','PDCD1','BTLA','TOX')
gene_SA_CD8_pre_effector = c('GZMK')
gene_SA_CD8_memory = c('ZNF683')
gene_SA_CD8_effector = c('NKG7', 'PRF1', 'CX3CR1')
gene_SA_CD8_Terminal_exhausted = c('CXCL13','PDCD1','BTLA','TOX')
gene_SA = c(gene_SA_CD4_memory,gene_SA_CD4_effector_memory,gene_SA_CD4_exhaustd,gene_SA_CD8_pre_effector,gene_SA_CD8_memory,gene_SA_CD8_effector,gene_SA_CD8_Terminal_exhausted)


##NM
library(readxl)
cellmarker_data <- read_excel("D:/nut/lab/TNBC_neo_scRNAseq/参考信息/cellmarker.xlsx")
gene_NM <- cellmarker_data[3:49,7:8]
colnames(gene_NM) <- c("GeneSymbol","class")


#######################
# ???????????????? FAnjia
expr_fanjia = as.data.frame(subcell_2@assays$RNA@data[gene_fanjia,])
# ????????????
meta = subcell_2@meta.data$seurat_clusters
names(meta)=rownames(subcell_2@meta.data)
meta = as.data.frame(meta)
colnames(meta)='subcluster'
# ??????cluster??????????????????????
new_expr_fanjia =do.call(rbind,
                         lapply(unique(meta$subcluster), function(t){
                           rowMeans(expr_fanjia[,meta$subcluster==t])
                         })
)
# ????????????
# ????????????
matrix_fanjia = as.matrix(new_expr_fanjia)
scale_matrix_fanjia = scale(matrix_fanjia) # scale??????????????????????????????????????????scale
scale_matrix_fanjia = t(scale_matrix_fanjia)
colnames(scale_matrix_fanjia) = unique(meta$subcluster)
# ??????marker?????????????? 
mylist_fanjia = list(gene_cell_type,gene_naive,gene_inhibitory,gene_cytotoxic,gene_co_stimulatory,gene_resident,gene_TF,gene_Treg)
names(mylist_fanjia) = c("cell","naive","inhibitory","cytotoxic","co_stimu","resident","TF","Treg")
mylist_fanjia_unlist =  as.data.frame (unlist(mylist_fanjia))  
colnames(mylist_fanjia_unlist)= 'gene'
mylist_fanjia_unlist$marker_type = rownames(mylist_fanjia_unlist)
mylist_fanjia_unlist$marker_type = str_sub(mylist_fanjia_unlist$marker_type,1,nchar(mylist_fanjia_unlist$marker_type)-1)
rownames(mylist_fanjia_unlist) <- mylist_fanjia_unlist$gene
# ????????????
library(ComplexHeatmap)
pdf(file = paste(sampleID,"Fanjia.pdf",sep = ""),height = 12,width = 8)
Heatmap(scale_matrix_fanjia,cluster_rows = F,split = mylist_fanjia_unlist$marker_type,cluster_columns = F )
dev.off()

pdf(file = paste(sampleID,"Fanjia_clu.pdf",sep = ""),height = 12,width = 8)
Heatmap(scale_matrix_fanjia,cluster_rows = F,split = mylist_fanjia_unlist$marker_type,cluster_columns = T )
dev.off()



#######################
# ???????????????? #gene_SA
expr_fanjia = as.data.frame(subcell_2@assays$RNA@data[gene_SA,])
# ????????????
meta = subcell_2@meta.data$seurat_clusters
names(meta)=rownames(subcell_2@meta.data)
meta = as.data.frame(meta)
colnames(meta)='subcluster'
# ??????cluster??????????????????????
new_expr_fanjia =do.call(rbind,
                         lapply(unique(meta$subcluster), function(t){
                           rowMeans(expr_fanjia[,meta$subcluster==t])
                         })
)
# ????????????
# ????????????
matrix_fanjia = as.matrix(new_expr_fanjia)
scale_matrix_fanjia = scale(matrix_fanjia) # scale??????????????????????????????????????????scale
scale_matrix_fanjia = t(scale_matrix_fanjia)
colnames(scale_matrix_fanjia) = unique(meta$subcluster)
# ??????marker?????????????? 
mylist_fanjia = list(gene_SA_CD4_memory,gene_SA_CD4_effector_memory,gene_SA_CD4_exhaustd,gene_SA_CD8_pre_effector,gene_SA_CD8_memory,gene_SA_CD8_effector,gene_SA_CD8_Terminal_exhausted)
names(mylist_fanjia) = c(" CD4_memory"," CD4_effector_memory"," CD4_exhaustd"," CD8_pre_effector",
                         " CD8_memory"," CD8_effector"," CD8_Terminal_exhausted")
mylist_fanjia_unlist =  as.data.frame (unlist(mylist_fanjia))  
colnames(mylist_fanjia_unlist)= 'gene'
mylist_fanjia_unlist$marker_type = rownames(mylist_fanjia_unlist)
mylist_fanjia_unlist$marker_type = str_sub(mylist_fanjia_unlist$marker_type,1,nchar(mylist_fanjia_unlist$marker_type)-1)

# ????????????
library(ComplexHeatmap)
pdf(file = paste(sampleID,"SA_gene_Anno_1.pdf",sep = ""),height = 12,width = 8)
Heatmap(scale_matrix_fanjia,cluster_rows = F,split = mylist_fanjia_unlist$marker_type,cluster_columns = F )
dev.off()



#######################
# ????????????????##NM????
expr_fanjia = as.matrix(subcell_2@assays$RNA@data[gene_NM$GeneSymbol,])
# ????????????
meta = subcell_2@meta.data$seurat_clusters
names(meta)=rownames(subcell_2@meta.data)
meta = as.data.frame(meta)
colnames(meta)='subcluster'
# ??????cluster??????????????????????
new_expr_fanjia =do.call(rbind,
                         lapply(unique(meta$subcluster), function(t){
                           rowMeans(expr_fanjia[,meta$subcluster==t])
                         })
)
# ????????????
# ????????????
matrix_fanjia = as.matrix(new_expr_fanjia)
scale_matrix_fanjia = scale(matrix_fanjia) # scale??????????????????????????????????????????scale
scale_matrix_fanjia = t(scale_matrix_fanjia)
colnames(scale_matrix_fanjia) = unique(meta$subcluster)

# ????????????
scale_matrix_fanjia_1 <- scale_matrix_fanjia[,as.character(c(0:12))]
library(ComplexHeatmap)
pdf(file = paste(sampleID,"NM_gene_Anno_1.pdf",sep = ""),height = 12,width = 8)
Heatmap(scale_matrix_fanjia,cluster_rows = F,split = gene_NM$class,cluster_columns =F,col = mycolors )
graphics.off()


pdf(file = paste(sampleID,"NM_gene_Anno_order.pdf",sep = ""),height = 12,width = 8)
Heatmap(scale_matrix_fanjia_1,cluster_rows = F,split = gene_NM$class,cluster_columns =F,col = mycolors )
graphics.off()

p <- Heatmap(scale_matrix_fanjia,cluster_rows = F,split = gene_NM$class,cluster_columns = T )
export::graph2pdf(p,file=paste(sampleID,"_NM_anno.pdf"),width = 8,height = 12)

########
##个别亚群看一看
tmp_marker <- FindMarkers(subcell_2,ident.1 = 14,ident.2 = 3)
write.csv(tmp_marker,file = "TNK_14vs3_findmarker.csv")

DimPlot(subcell_2,cells = 14)
FeaturePlot(subcell_2,features = c("CXCL8", "FCGR3B", "MNDA", "SELL"),cells = WhichCells(subcell_2,idents = 15))
FeaturePlot(subcell_2,features = c("MALAT1"),cells = WhichCells(subcell_2,idents = 15))
FeaturePlot(subcell_2,features = c("IGLC2", "IGKC"),cells = WhichCells(subcell_2,idents = 2))
VlnPlot(subcell_2,features = c("CD3E","IGHG1","MALAT1"),pt.size = 0)
VlnPlot(subcell_2,features = c("CD3E","IGHG1","NEAT1","percent.mt"),pt.size = 0)

colnames(subcell_2@meta.data)
########################
###????????

subcell_2 <- RenameIdents(subcell_2,
                          `0` = "CD8T_GZMK",
                          `1` = "CD4T_IL7R",
                          `2` = "Unassigned1",
                          `3` = "CD4T_FOXP3",
                          `4` = "CD8T_GZMB",
                          `5` = "T_IFIT3",
                          `6` = "gdT_TRDV2",
                          `7` = "CD4_CXCL13",
                          `8` = "T_CCR7",
                          `9` = "NK_XCL1",
                          `10` = "NK_FGFBP2",
                          `11` = "T_MKI67",
                          `12` = "Unassigned2")


subcell_2 <- AddMetaData(subcell_2,subcell_2@active.ident,col.name = "T_cell_recluster_anno")



subcell_2@meta.data$Response_all[subcell_2$orig.ident_anno%in%pCR] <- "pCR"
subcell_2@meta.data$Response_all[subcell_2$orig.ident_anno%in%non_PCR] <- "non_pCR"

subcell_2@meta.data$Response_detail[subcell_2$orig.ident_anno%in%Immu_pCR] <- "ImmupCR"
subcell_2@meta.data$Response_detail[subcell_2$orig.ident_anno%in%Immu_nonpCR] <- "Immu_nonpCR"
subcell_2@meta.data$Response_detail[subcell_2$orig.ident_anno%in%Chem_pCR] <- "Chem_pCR"
subcell_2@meta.data$Response_detail[subcell_2$orig.ident_anno%in%Chem_nonpCR] <- "Chem_nonpCR"



save(subcell_2,file = "T_cell_anno.Rdata")
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
col_T_16 = getPalette(13)

p <- DimPlot(subcell_2,cols = col_T_16,raster = F)
p
ggsave(plot = p,filename = paste(sampleID,"_dimplot.pdf"),width = 8,height = 6)
#export::graph2ppt(p,file=paste(sampleID,"_dimplot.pptx"),width = 8,height = 6)
load("D:/nut/lab/TNBC_neo_scRNAseq/36sample_final/Step2.all_cell_anno/V2/sample_order.Rdata")
GroupID <- names(table(subcell_2$orig.ident_anno))
result <- data.frame(
  celltype = names(table(subcell_2@meta.data$T_cell_recluster_anno))
)
for (i in 1:length(GroupID)) {
  Data <- subcell_2@meta.data[subcell_2@meta.data$orig.ident_anno==GroupID[i],]
  freq <- table(Data$T_cell_recluster_anno)
  percent = as.data.frame(round(freq/sum(freq)*100, 2))
  colnames(percent) <- c("celltype",GroupID[i])
  result <- left_join(result,percent,by="celltype")
}
df <- subcell_2@meta.data
df$orig.ident_anno <- factor(df$orig.ident_anno,levels = sample_order )
p <- df %>% ggplot(aes(x = orig.ident_anno, fill = T_cell_recluster_anno)) + #x??????????clarity????????????color??J??H?? 
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =col_T_16) + #?????????? 
  theme_classic() + #???????? 
  labs(y = 'Percent') 
p
ggsave(plot = p,filename = paste(sampleID,"_fraction_patient.pdf"),width = 12,height = 8)
#export::graph2ppt(p,file=paste(sampleID,"_raction_patient.pptx"),width = 6,height = 8)


P <-df %>% ggplot(aes(x = Time, fill = T_cell_recluster_anno)) + #x轴的分类为clarity，填充颜色为color（J和H） 
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =col_T_16) + #设置颜色板 
  theme_classic() + #设置主题 
  labs(y = 'Percent') +
  labs(x='')
#+ #设置y轴名为‘Percent’ 
ggsave(plot = P,filename = paste(sampleID,'_prop_of_T_cell_recluster_anno_in_Time.pdf'),width = 6,height = 8)

P <-df %>% ggplot(aes(x = Stage, fill = T_cell_recluster_anno)) + #x轴的分类为clarity，填充颜色为color（J和H） 
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =col_T_16) + #设置颜色板 
  theme_classic() + #设置主题 
  labs(y = 'Percent') +
  labs(x='')
#+ #设置y轴名为‘Percent’ 
ggsave(plot = P,filename = paste(sampleID,'_prop_of_T_cell_recluster_anno_in_Stage.pdf'),width = 6,height = 8)

P <-df %>% ggplot(aes(x = Response_detail, fill = T_cell_recluster_anno)) + #x轴的分类为clarity，填充颜色为color（J和H） 
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =col_T_16) + #设置颜色板 
  theme_classic() + #设置主题 
  labs(y = 'Percent') +
  labs(x='')

P#+ #设置y轴名为‘Percent’ 
ggsave(plot = P,filename = paste(sampleID,'_prop_of_T_cell_recluster_anno_in_response.pdf'),width = 12,height = 8)

P <-df %>% ggplot(aes(x = Response_all, fill = T_cell_recluster_anno)) + #x轴的分类为clarity，填充颜色为color（J和H） 
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =col_T_16) + #设置颜色板 
  theme_classic() + #设置主题 
  labs(y = 'Percent') +
  labs(x='')

P#+ #设置y轴名为‘Percent’ 
ggsave(plot = P,filename = paste(sampleID,'_prop_of_T_cell_recluster_anno_in_chem_response.pdf'),width = 12,height = 8)



###########临床关注的重要免疫调节
library(Seurat)
library(MAST)
library(insight)
library(ggpubr)
library(ggplot2)
library(future)
library(future.apply)
library(stringr)

options(future.globals.maxSize = 25000 * 1048^2)
plan("multisession", workers = 8)
options(scipen=10000)

map_colours <- c("#b2182b", "#1F618D",
                 "#f4a582", "#9B59B6","#85929E", "#1c9099", "#74add1", "#053061",
                 "#1b7837", "#b8e186", "#bebada", "#fed976",
                 "#e7298a", "#47E3FF", "#FFBF47",
                 "#A93226", "#270e26","orange",
                 "#b8bc53", "#5628ce", "#fa909c",
                 "#8ff331", "#FF6347", "#6347FF", "#556270", "#4ECDC4", "#C7F464", "#FF6B6B", "#C44D58", "#E3FF47","#FF4787")


#Read files
genes_list_file <- read.csv("D:/nut/lab/TNBC_neo_scRNAseq/参考信息/myloid_marker/clinical_immunegene_of_interest.csv")


#Subset out Myeloid and T-cells, split by celltype then run DGE for each BrCa subtype comparison
DefaultAssay(subcell_2) <- "RNA"
Idents(object = subcell_2) <- "seurat_clusters"
sample_list <- SplitObject(subcell_2, split.by = "seurat_clusters")
table(subcell_2$seurat_clusters)
aaa <- data.frame()
for(i in 1:length(sample_list)){
  tryCatch({
    name <- names(sample_list[i])
    aaa <- FindMarkers(sample_list[[i]], group.by = "Response_detail",ident.1 = "ImmupCR",ident.2 = "Immu_nonpCR", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    bbb <- FindMarkers(sample_list[[i]], group.by = "Response_detail",ident.1 = "Chem_pCR",ident.2 = "Chem_nonpCR", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    ccc <- FindMarkers(sample_list[[i]], group.by = "Response_all",ident.1 = "pCR",ident.2 = "non_pCR", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    ddd <- FindMarkers(sample_list[[i]], group.by = "Time",ident.1 = "Pre_neo",ident.2 = "Advanced", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    write.csv(aaa, file = paste0(name,"_MAST_across_Immune_pCR_nonpCR.csv"), row.names = TRUE)
    write.csv(bbb, file = paste0(name,"_MAST_across_Chem_pCR_nonpCR.csv"), row.names = TRUE)
    write.csv(ccc, file = paste0(name,"_MAST_across_pCR_nonpCR.csv"), row.names = TRUE)
    write.csv(ddd, file = paste0(name,"_MAST_across_Stage_early_adv.csv"), row.names = TRUE)
  }, error=function(e){print("error")})
}


#For TNBC vs HER2 comparison, generate dataframe of the genes of interest for each celltype comparison
path_file <- "./DEG/"
dge_files <- list.files(path = path_file)
gene_OI <- genes_list_file$gene_name
gene_OI <- c(gene_OI,"MERTK","AXL","TYRO3")
file_check <- read.csv(file = paste0(path_file,dge_files[1]),row.names = 1)
checkpoints_significant <- data.frame(matrix(ncol = 5, nrow = 0))


for (i in 1:length(dge_files)){
  tryCatch({
    file_check <- read.csv(file = paste0(path_file,dge_files[i]),row.names = 1)
    name_to_use <- gsub("_MAST_across", "", dge_files[i])
    name_to_use <- gsub(".csv", "", name_to_use)
    tmp <- unlist(str_split(name_to_use,fixed("_")))
    inter_df <- file_check[which(rownames(file_check) %in% gene_OI),]
    inter_df$celltype <- tmp[1]
    inter_df$comparison <- tmp[2]
    inter_df$gene <- paste(tmp[length(tmp)-2],tmp[length(tmp)-1],tmp[length(tmp)],rownames(inter_df),sep="_")
    checkpoints_significant <- rbind(checkpoints_significant,inter_df)
  }, error=function(e){print("error")})
}


C2 <- checkpoints_significant


#Aggregate all data.frames into one and add NS for any none signficant values
all_data <- checkpoints_significant
all_data$sign <- insight::format_p(all_data$p_val_adj,stars_only = TRUE)
all_data <- all_data %>% mutate(sign = ifelse(all_data$p_val_adj>0.05,"ns",sign))
write.csv(all_data, file = "paired_comparison_DGE_MAST.csv", row.names = T)

#Arrange by high to low Log threshold
plot_data <- all_data
plot_data$gene <- as.factor(plot_data$gene)
plot_data$comparison <- as.factor(plot_data$comparison)
plot_data <- plot_data %>% group_by(comparison) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
plot_data <- plot_data %>% arrange(desc(avg_log2FC))

#Remove all not signicant
plot_data2 <- subset(plot_data,!plot_data$sign=="ns")


#best figure, but other options are out there
a <- ggdotchart(plot_data2, x = "gene", y = "avg_log2FC",
                color = "celltype",                          # Color by groups
                palette = map_colours, # Custom color palette
                sorting = "descending",
                add = "segments",                # Sort value in descending order                    # Add segments from y = 0 to dots
                add.params = list(color = "lightgray", size = 1), # Change segment color and size
                group = "comparison",                                # Order by groups
                dot.size = 6,                                 # Large dot size
                # label = round(All_signf$p_val_adj,1),                        # Add mpg values as dot labels
                label = "sign",             # Add mpg values as dot labels
                font.label = list(color = "white", size = 7,vjust = 0.5),               # Adjust label parameters
                ggtheme = theme_pubr()) +
  theme_cleveland() + geom_hline(yintercept = 0, linetype = 1, color = "black") +
  coord_cartesian(ylim = c(-1, 1)) + facet_wrap(~comparison, nrow=1, scales="free_x")

a <- a+font("xy.text", size = 9, color = "black", face = "bold.italic")
a
ggsave(plot = a,filename = "myeloid_ICB_lollipo.pdf",width = 12,height = 6)

#####-----immune_modulator
library(Seurat)
library(tidyverse)
library(viridis)
library(pheatmap)

genes_OI <- read.csv("X:scRNAseq_refData/2021-NG/cODEBrCa_cell_atlas-main/downstream_immune_analysis/gene_list/immune_modulator_gene_list.csv")
genes <- intersect(subcell_2@assays[["RNA"]]@data@Dimnames[[1]],as.character(unique(genes_OI$锘縢ene_name)))
dplot <- DotPlot(object = subcell_2, features = genes, split.by = "Response_detail",cols = rainbow(12))


subcell_3 <- subset(subcell_2,Immune_res=="Good"|Immune_res=="Bad")
dplot <- DotPlot(object = subcell_3, features = genes, split.by = "Immune_res")


ddata <- as.data.frame(dplot[["data"]])
df_spread <- tidyr::spread(ddata[, c(3,4,5)], id, avg.exp.scaled)
rownames(df_spread) <- df_spread[, 1]
df_spread <- df_spread[, -1]
mat <- as.matrix(df_spread)
mat <- na.omit(mat)

phet<- pheatmap(mat,
                color = inferno(10),
                fontsize = 6,
                cellheight = 6,
                border_color = T,
                fontsize_row = 6,
                cellwidth = 6)

phet
ggsave(phet, filename = "pheatmap_all_BrCa_NKT.pdf", width = 8, height = 16)



fucnction_score <- read.csv("X:/scRNAseq_refData/2021-NG/cODEBrCa_cell_atlas-main/downstream_immune_analysis/gene_list/fucnction_score.csv")
colnames(fucnction_score)[1] <- "Dysfunction_score"
Dysfunction_score <- list(c(fucnction_score$Dysfunction_score))
Cytotoxic_score <- intersect(fucnction_score$Cytotoxic_score,rownames(subcell_Mono@assays[["RNA"]]@counts))
Cytotoxic_score <- list(c(Cytotoxic_score))
Cytoyoxic_gene_NG <- intersect(fucnction_score$Cytoyoxic_gene_NG,rownames(subcell_Mono@assays[["RNA"]]@counts))
Cytoyoxic_gene_NG <- list(c(Cytoyoxic_gene_NG))

subcell_Mono <-AddModuleScore( object = subcell_Mono, features = Dysfunction_score, ctrl = 100, name = 'Dysfunction_score' )
subcell_Mono <-AddModuleScore( object = subcell_Mono, features = Cytotoxic_score, ctrl = 100, name = 'Cytotoxic_score' )
subcell_Mono <-AddModuleScore( object = subcell_Mono, features = Cytoyoxic_gene_NG, ctrl = 100, name = 'Cytoyoxic_gene_NG')
colnames(subcell_Mono@meta.data)
colnames(subcell_Mono@meta.data)[43:45] <- c('Dysfunction_score','Cytotoxic_score','Cytoyoxic_gene_NG' ) 


p <- VlnPlot(subcell_Mono,features = c('Dysfunction_score'),pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
p
ggsave(plot = p,filename = "T_Dysfunction_score_vln.pdf",width = 8,height = 4)

p <- VlnPlot(subcell_Mono,features = c('Cytoyoxic_gene_NG'),pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ NoLegend()
p
ggsave(plot = p,filename = "T_Cytoyoxic_score_vln.pdf",width = 8,height = 4)

table(subcell_Mono$T_cell_recluster_anno)
p <- DimPlot(subcell_Mono,cols = col_T_16,group.by = "T_cell_recluster_anno",raster = T)
p
ggsave(plot = p,filename = paste(sampleID,"_dimplot.pdf"),width = 8,height = 6)
