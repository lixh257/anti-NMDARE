 FigS1

1. Fig S1A(Mitosort),展示分的样本和预测的doublets,以及其中计算得到的P值
#source:MitoSort.sh
##heatmap &p value plot


setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
###Lib1
#$cp -r /data/R03/zhangwx/project/separeteLib1/MitoSort/ /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib1/
pre_Lib1=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib1/MitoSort/Demultiplex_output/result_pvalue.txt",sep="\t",header=T)
pre_Lib1_patient1=gsub("-",".",pre_Lib1[which(pre_Lib1[,2]=="Sample0"),1])
pre_Lib1_patient2=gsub("-",".",pre_Lib1[which(pre_Lib1[,2]=="Sample1"),1])
pre_Lib1_doublet=gsub("-",".",pre_Lib1[which(pre_Lib1[,2]=="Doublet"),1])
barcode_order=c(pre_Lib1_patient1,pre_Lib1_patient2,pre_Lib1_doublet)
mat=read.csv("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib1/MitoSort/SNP_matrix/frequency.csv",head=T,row.names=1,fill=T)
germline_SNP=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib1/MitoSort/Demultiplex_output/specific_germline.txt",sep="\t",head=T)
germline_SNP_n=germline_SNP[,1]
mat[is.na(mat)]=0
mat1=mat[which(rownames(mat) %in% germline_SNP_n),barcode_order]
##按照SNP是哪个样本的计算mat中mutation rate 的均值
snp1=germline_SNP[which(germline_SNP[,2]=="Sample0"),1]
snp2=germline_SNP[which(germline_SNP[,2]=="Sample1"),1]
#mat1=as.matrix(mat1)
mat1_mean_1=apply(mat1[which(rownames(mat1)%in%snp1),],2,mean)
mat1_mean_2=apply(mat1[which(rownames(mat1)%in%snp2),],2,mean)
mat_res=rbind(mat1_mean_1,mat1_mean_2)
rownames(mat_res)=c("Patient1_SNP","Patient2_SNP")
write.table(mat_res,"Mitosort_res_Lib1_matrix_mean.txt",sep="\t",row.names=T,col.names=T,quote=F)

mat=read.table("Mitosort_res_Lib1_matrix_mean.txt",sep="\t",head=T,row.names=1)
Sample_color <- paletteer::paletteer_d("ggsci::category20_d3")
ann_col=data.frame(Sample=factor(c(rep("Patient1",length(pre_Lib1_patient1)),
                                   rep("Patient2",length(pre_Lib1_patient2)),
                                   rep("Doublet",length(pre_Lib1_doublet)))))
annoColor <- list(Sample=c("Patient1"=Sample_color[11],"Patient2"=Sample_color[12],"Doublet"=Sample_color[15]))
colnames(mat)=rownames(ann_col)
library(pheatmap)
pdf("FigS1B_Lib1_mitosort_res.pdf",width=8,height=3)
pheatmap(mat,    
color = colorRampPalette(paletteer::paletteer_d("PNWColors::Shuksan2"))(10),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = FALSE,   
annotation_col = ann_col,  annotation_colors = annoColor
)
dev.off()

###########################################################################################
#Lib3_add:
#com=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter.rds")
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
mat=read.csv("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib3/MitoSort/SNP_matrix/frequency.csv",head=T,row.names=1,fill=T)
germline_SNP=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib3/MitoSort/Demultiplex_output/specific_germline.txt",sep="\t",head=T)
pre_Lib3=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib3/MitoSort/Demultiplex_output/lib3_depth_1_result_pvalue.txt",sep="\t",head=T)
#pre_res[,1]=gsub("-",".",pre_res[,1])
pre_Lib3_patient=gsub("-",".",pre_Lib3[which(pre_Lib3[,2]=="Sample1"),1])
pre_Lib3_Normal1=gsub("-",".",pre_Lib3[which(pre_Lib3[,2]=="Sample0"),1])
pre_Lib3_doublet=gsub("-",".",pre_Lib3[which(pre_Lib3[,2]=="Doublet"),1])
pre_Lib3_unassign=gsub("-",".",pre_Lib3[which(pre_Lib3[,2]=="Unassign"),1])
barcode_order=c(pre_Lib3_patient,pre_Lib3_Normal1,pre_Lib3_unassign,pre_Lib3_doublet)
germline_SNP_n=germline_SNP[,1]
mat[is.na(mat)]=0
mat1=mat[which(rownames(mat) %in% germline_SNP_n),barcode_order]
##按照SNP是哪个样本的计算mat中mutation rate 的均值

snp1=germline_SNP[which(germline_SNP[,2]=="Sample0"),1]
snp2=germline_SNP[which(germline_SNP[,2]=="Sample1"),1]
#mat1=as.matrix(mat1)
mat1_mean_1=apply(mat1[which(rownames(mat1)%in%snp1),],2,mean)
mat1_mean_2=apply(mat1[which(rownames(mat1)%in%snp2),],2,mean)
mat_res=rbind(mat1_mean_1,mat1_mean_2)
rownames(mat_res)=c("Patient3_SNP","Normal1_SNP")
write.table(mat_res,"Mitosort_res_Lib3_matrix_mean.txt",sep="\t",row.names=T,col.names=T,quote=F)

mat=read.table("Mitosort_res_Lib3_matrix_mean.txt",sep="\t",head=T,row.names=1)
Sample_color <- paletteer::paletteer_d("ggsci::category20_d3")
ann_col=data.frame(Sample=factor(c(rep("Patient3",length(pre_Lib3_patient)),
                                   rep("Normal1",length(pre_Lib3_Normal1)),
                                   rep("Unassign",length(pre_Lib3_unassign)),
                                   rep("Doublet",length(pre_Lib3_doublet))
                                   )))
annoColor <- list(Sample=c("Patient3"=Sample_color[11],"Normal1"=Sample_color[12],
                           "Doublet"=Sample_color[15],"Unassign"=Sample_color[9]))
colnames(mat)=rownames(ann_col)
library(pheatmap)
pdf("FigS1B_Lib3_mitosort_res.pdf",width=8,height=3)
pheatmap(mat,    
color = colorRampPalette(paletteer::paletteer_d("PNWColors::Shuksan2"))(10),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = FALSE,   
annotation_col = ann_col,  annotation_colors = annoColor
)
dev.off()


############################################################################################
#Lib4_add:
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
mat=read.csv("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib4/MitoSort/SNP_matrix/frequency.csv",head=T,row.names=1,fill=T)
germline_SNP=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib4/MitoSort/Demultiplex_output/specific_germline.txt",sep="\t",head=T)
pre_Lib4=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib4/MitoSort/Demultiplex_output/lib4_depth_1_result_pvalue.txt",sep="\t",head=T)

#pre_res[,1]=gsub("-",".",pre_res[,1])
pre_Lib4_patient4=gsub("-",".",pre_Lib4[which(pre_Lib4[,2]=="Sample1"),1])
pre_Lib4_Normal2=gsub("-",".",pre_Lib4[which(pre_Lib4[,2]=="Sample0"),1])
pre_Lib4_doublet=gsub("-",".",pre_Lib4[which(pre_Lib4[,2]=="Doublet"),1])
pre_Lib4_unassign=gsub("-",".",pre_Lib4[which(pre_Lib4[,2]=="Unassign"),1])
barcode_order=c(pre_Lib4_patient4,pre_Lib4_Normal2,pre_Lib4_unassign,pre_Lib4_doublet)
germline_SNP_n=germline_SNP[,1]
mat[is.na(mat)]=0
mat1=mat[which(rownames(mat) %in% germline_SNP_n),barcode_order]
##按照SNP是哪个样本的计算mat中mutation rate 的均值
snp1=germline_SNP[which(germline_SNP[,2]=="Sample0"),1]
snp2=germline_SNP[which(germline_SNP[,2]=="Sample1"),1]
#mat1=as.matrix(mat1)
mat1_mean_1=apply(mat1[which(rownames(mat1)%in%snp1),],2,mean)
mat1_mean_2=apply(mat1[which(rownames(mat1)%in%snp2),],2,mean)
mat_res=rbind(mat1_mean_1,mat1_mean_2)
rownames(mat_res)=c("Patient4_SNP","Normal2_SNP")
write.table(mat_res,"Mitosort_res_Lib4_matrix_mean.txt",sep="\t",row.names=T,col.names=T,quote=F)

mat=read.table("Mitosort_res_Lib4_matrix_mean.txt",sep="\t",head=T,row.names=1)
Sample_color <- paletteer::paletteer_d("ggsci::category20_d3")
ann_col=data.frame(Sample=factor(c(rep("Patient4",length(pre_Lib4_patient4)),
                                   rep("Normal2",length(pre_Lib4_Normal2)),
                                   rep("Unassign",length(pre_Lib4_unassign)),
                                   rep("Doublet",length(pre_Lib4_doublet))
                                   )))
annoColor <- list(Sample=c("Patient4"=Sample_color[11],"Normal2"=Sample_color[12],
                           "Doublet"=Sample_color[15],"Unassign"=Sample_color[9]))
colnames(mat)=rownames(ann_col)
library(pheatmap)
pdf("FigS1B_Lib4_mitosort_res.pdf",width=8,height=3)
pheatmap(mat,    
color = colorRampPalette(paletteer::paletteer_d("PNWColors::Shuksan2"))(10),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = FALSE,   
annotation_col = ann_col,  annotation_colors = annoColor
)
dev.off()


###Lib5
##cp -r /data/R01/Chenzh275/Project/Zhangwx/lib5/ /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib5/
mat=read.csv("/data/R01/Chenzh275/Project/Zhangwx/lib5/MitoSort/SNP_matrix/frequency.csv",head=T,row.names=1,fill=T)
germline_SNP=read.table("/data/R01/Chenzh275/Project/Zhangwx/lib5/MitoSort/Demultiplex_output/specific_germline.txt",sep="\t",head=T)
pre_Lib5=read.table("/data/R01/Chenzh275/Project/Zhangwx/lib5/MitoSort/Demultiplex_output/result_pvalue.txt",sep="\t",header=T)
pre_Lib5_normal3=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Sample0"),1])
pre_Lib5_normal4=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Sample1"),1])
pre_Lib5_normal5=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Sample2"),1])
pre_Lib5_normal6=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Sample3"),1])
pre_Lib5_doublet=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Doublet"),1])
pre_Lib5_unassign=gsub("-",".",pre_Lib5[which(pre_Lib5[,2]=="Unassign"),1])
barcode_order=c(pre_Lib5_normal3,pre_Lib5_normal4,pre_Lib5_normal5,pre_Lib5_normal6,pre_Lib5_unassign,pre_Lib5_doublet)
germline_SNP_n=germline_SNP[,1]
mat[is.na(mat)]=0
mat1=mat[which(rownames(mat) %in% germline_SNP_n),barcode_order]
##按照SNP是哪个样本的计算mat中mutation rate 的均值
###这个specific_germline.txt里面的ID和mat里面的不符合
snp1=germline_SNP[which(germline_SNP[,2]=="Sample1"),1]
snp2=germline_SNP[which(germline_SNP[,2]=="Sample2"),1]
snp3=germline_SNP[which(germline_SNP[,2]=="Sample0"),1]
snp4=germline_SNP[which(germline_SNP[,2]=="Sample3"),1]
#mat1=as.matrix(mat1)
mat1_mean_1=apply(mat1[which(rownames(mat1)%in%snp1),],2,mean)
mat1_mean_2=apply(mat1[which(rownames(mat1)%in%snp2),],2,mean)
mat1_mean_3=apply(mat1[which(rownames(mat1)%in%snp3),],2,mean)
mat1_mean_4=apply(mat1[which(rownames(mat1)%in%snp4),],2,mean)
mat_res=rbind(mat1_mean_1,mat1_mean_2,mat1_mean_3,mat1_mean_4)
rownames(mat_res)=c("Normal3_SNP","Normal4_SNP","Normal5_SNP","Normal6_SNP")
write.table(mat_res,"Mitosort_res_Lib5_matrix_mean.txt",sep="\t",row.names=T,col.names=T,quote=F)

mat=read.table("Mitosort_res_Lib5_matrix_mean.txt",sep="\t",head=T,row.names=1)
Sample_color <- paletteer::paletteer_d("ggsci::category20_d3")
ann_col=data.frame(Sample=factor(c(rep("Normal3",length(pre_Lib5_normal3)),
                                   rep("Normal4",length(pre_Lib5_normal4)),
                                   rep("Normal5",length(pre_Lib5_normal5)),
                                   rep("Normal6",length(pre_Lib5_normal6)),
                                   rep("Unassign",length(pre_Lib5_unassign)),
                                   rep("Doublet",length(pre_Lib5_doublet))
                                   )))

annoColor <- list(Sample=c("Normal3"=Sample_color[11],"Normal4"=Sample_color[12],
                           "Normal5"=Sample_color[13],"Normal6"=Sample_color[14],
                           "Doublet"=Sample_color[15],"Unassign"=Sample_color[9]))                        
colnames(mat)=rownames(ann_col)
library(pheatmap)
pdf("FigS1B_Lib5_mitosort_res.pdf",width=8,height=4)
pheatmap(mat,    
color = colorRampPalette(paletteer::paletteer_d("PNWColors::Shuksan2"))(10),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = FALSE,   
annotation_col = ann_col,  annotation_colors = annoColor
)
dev.off()


1. Fig S1B(QC plot)
#BiocManager::install("dittoSeq")
library(dittoSeq)
library(dplyr)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
Lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
Lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")
Lib3=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")
Lib4=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")
Lib5=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")

Lib1$dataset <- 'Lib1' 
Lib2$dataset <- 'Lib2' 
Lib3$dataset <- 'Lib3' 
Lib4$dataset <- 'Lib4' 
Lib5$dataset <- 'Lib5'
Merge_befor_filter <- merge(   x = Lib1,   y = list(Lib2,Lib3,Lib4,Lib5),   add.cell.ids = c("Lib1","Lib2", "Lib3", "Lib4","Lib5") )
saveRDS(Merge_befor_filter,"Merge_befor_filter.rds")


############################################################
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
Merge_befor_filter=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_befor_filter.rds")
1.6.1  Fig S1A(QC_vlnplot)
#Fraction of reads in peaks (FRiP) 
pdf("QC_vlnplot_boxplot_FRiP.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "FRiP", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'the fraction of reads in peaks',
          y.breaks = seq(0,1,0.2),
          xlab = '',
         # x.labels = c("Female","Male"),
          x.labels.rotate =F,
          max=1,
          min=0,
          main = "FRiP",
          legend.show = F)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
pdf("QC_vlnplot_boxplot_FRiP_Lib4.pdf",width=3,height=5)
dittoPlot(Lib4,
          "FRiP", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'the fraction of reads in peaks',
          y.breaks = seq(0,1,0.2),
          xlab = '',
         # x.labels = c("Female","Male"),
          x.labels.rotate =F,
          max=1,
          min=0,
          main = "FRiP",
          legend.show = F)
dev.off()

#MT
pdf("QC_vlnplot_boxplot_MTreads.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "percent.mt", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'MT reads precentage',
          y.breaks = seq(0,20,1),
          xlab = '',
          x.labels.rotate =F,
          max=20,
          min=0,
          main = "MT reads",
          legend.show = F)
dev.off()

#nucleosome_signal
pdf("QC_vlnplot_boxplot_nucleosome_signal.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "nucleosome_signal", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'nucleosome signal',
          y.breaks = seq(0,4,1),
          xlab = '',
          x.labels.rotate =F,
          max=4,
          min=0,
          main = "nucleosome_signal",
          legend.show = F)
dev.off()

#nCount_RNA&'nCount_ATAC
pdf("QC_vlnplot_boxplot_nCount_RNA.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "nCount_RNA", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'nCount_RNA',
          y.breaks = seq(0,10000,500),
          xlab = '',
          x.labels.rotate =F,
          max=10000,
          min=0,
          main = "nCount_RNA",
          legend.show = F)
dev.off()
pdf("QC_vlnplot_boxplot_nCount_ATAC.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "nCount_ATAC", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'nCount_ATAC',
          y.breaks = seq(0,50000,2000),
          xlab = '',
          x.labels.rotate =F,
          max=50000,
          min=0,
          main = "nCount_ATAC",
          legend.show = F)
dev.off()

#TSS enrichment
pdf("QC_vlnplot_boxplot_TSS.enrichment.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "TSS.enrichment", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'TSS.enrichment',
          y.breaks = seq(0,20,1),
          xlab = '',
          x.labels.rotate =F,
          max=20,
          min=0,
          main = "TSS",
          legend.show = F)
dev.off()

#others
#atac_fragments
pdf("QC_vlnplot_boxplot_atac_fragments.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "atac_fragments", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'atac fragments',
          y.breaks = seq(0,60000,10000),
          xlab = '',
          x.labels.rotate =F,
          max=60000,
          min=0,
          main = "atac_fragments",
          legend.show = F)
dev.off()

#UMI(gex_umis_count)
pdf("QC_vlnplot_boxplot_gex_umis_count.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "gex_umis_count", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'UMI counts',
          y.breaks = seq(0,15000,5000),
          xlab = '',
          x.labels.rotate =F,
          max=15000,
          min=0,
          main = "gex_umis_count",
          legend.show = F)
dev.off()

#gex_genes
pdf("QC_vlnplot_boxplot_gex_genes_count.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "gex_genes_count", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'Gene counts',
          y.breaks = seq(0,6000,1000),
          xlab = '',
          x.labels.rotate =F,
          max=6000,
          min=0,
          main = "gex_genes_count",
          legend.show = F)
dev.off()
#gex/atac_raw_reads(seq-depth)
pdf("QC_vlnplot_boxplot_atac_raw_reads.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "atac_raw_reads", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'atac_raw_reads',
          y.breaks = seq(0,200000,20000),
          xlab = '',
          x.labels.rotate =F,
          max=200000,
          min=0,
          main = "atac_raw_reads",
          legend.show = F)
dev.off()
pdf("QC_vlnplot_boxplot_gex_raw_reads.pdf",width=6,height=5)
dittoPlot(Merge_befor_filter ,
          "gex_raw_reads", 
          group.by = "dataset",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3,5,7,9),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 0.5),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'gex_raw_reads',
          y.breaks = seq(0,200000,20000),
          xlab = '',
          x.labels.rotate =F,
          max=200000,
          min=0,
          main = "gex_raw_reads",
          legend.show = F)
dev.off()




3. Fig S1C
###FigS1C cellnumber barplot
library(Seurat)
library(ggplot2)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
data=table(combined$dataset)
data=as.data.frame(data)
colnames(data)=c("Lib","Cellnumber")

pdf("FigS1C_cellnumber_Library.pdf",width=8,height=6)
ggplot(data,
aes(x=factor(Lib,levels=c("Lib1","Lib2","Lib3","Lib4","Lib5")),
y=Cellnumber,fill=Lib))+
geom_bar(position="stack",stat="identity")+
theme_bw()+
labs(x="library",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()



4.Fig S1D (cell distribution)

##cell ratio 
pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#539092","#8BCFCC","#51ADCF","#0278AE","#01A9B4","#A3D8F4")
colnames(combined@meta.data)
cellnumber=apply(table(combined@meta.data[,c(48,60)]),1,sum)
data=table(combined@meta.data[,c(48,60)])/apply(table(combined@meta.data[,c(48,60)]),1,sum)
data1=table(combined@meta.data[,60])/length(combined@meta.data[,60])
data_P=table(combined@meta.data[which(combined$AE=="Patient"),60])/length(combined@meta.data[which(combined$AE=="Patient"),60])
data_N=table(combined@meta.data[which(combined$AE=="Normal"),60])/length(combined@meta.data[which(combined$AE=="Normal"),60])
data_P=as.data.frame(data_P)
data_N=as.data.frame(data_N)
data=as.data.frame(data)
data1=as.data.frame(data1)

data1=cbind(rep("All_samples",length(data1[,1])),data1)
colnames(data1)=c("Label","celltype2","Freq")
data_P=cbind(rep("Patient",length(data_P[,1])),data_P)
colnames(data_P)=c("Label","celltype2","Freq")
data_N=cbind(rep("Normal",length(data_N[,1])),data_N)
colnames(data_N)=c("Label","celltype2","Freq")

data2=rbind(data,data1,data_P,data_N)
pal <- paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20, 3,13,19, 7,15,17,5,14,11, 12,16,18,8)]
pdf("Fig1D_barplot_cellratio.pdf",width=8,height=6)
ggplot(data2,
aes(x=factor(Label,levels=c("Patient1","Patient3","Patient4","Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","All_samples","Patient","Normal")),
y=Freq,
fill=factor(celltype2,levels=c("B naive","B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK","CM","ncM","cDC","pDC"))))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="proportion")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()



###Cell number
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
colnames(combined@meta.data)
cell_n=table(combined@meta.data[,c(48,60)])
cell_n=as.data.frame(cell_n)
colnames(cell_n)=c("Label","celltype2","Cellnumber")
data_P=table(combined@meta.data[which(combined$AE=="Patient"),60])
data_N=table(combined@meta.data[which(combined$AE=="Normal"),60])
data_P=as.data.frame(data_P)
data_N=as.data.frame(data_N)

data1=table(combined@meta.data[,60])
data1=as.data.frame(data1)
data1=cbind(rep("All_samples",length(data1[,1])),data1)
colnames(data1)=c("Label","celltype2","Cellnumber")
data_P=cbind(rep("Patient",length(data_P[,1])),data_P)
colnames(data_P)=c("Label","celltype2","Cellnumber")
data_N=cbind(rep("Normal",length(data_N[,1])),data_N)
colnames(data_N)=c("Label","celltype2","Cellnumber")

data2=cell_n
pal <- paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20, 3,13,19, 7,15,17,5,14,11, 12,16,18,8)]
pdf("FigS1D_barplot_cellnumber.pdf",width=8,height=6)
ggplot(data2,
aes(x=factor(Label,levels=c("Patient1","Patient3","Patient4","Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","All_samples","Patient","Normal")),
y=Cellnumber,
fill=factor(celltype2,levels=c("B naive","B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK","CM","ncM","cDC","pDC"))))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()



5.Fig S1E
(按照每个样本中cell的比例计算Bray-Curtis distance)
library(Seurat)
library(ggplot2)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
colnames(combined@meta.data)
data=table(combined@meta.data[,c(48,60)])/apply(table(combined@meta.data[,c(48,60)]),2,sum)
res=matrix(0,ncol=length(data[,1]),nrow=length(data[,1]))
colnames(res)=rownames(data)
rownames(res)=rownames(data)
data1=c()
for(i in 1:length(data[,1])){
for(j in 1:length(data[,1])){
data1=rbind(data[i,],data[j,])
res[i,j]=1-sum(apply(data1, 2, function(x) abs(max(x)-min(x)))) / sum(rowSums(data1))
}
}
library(pheatmap)
pheatmap(res,scale = "none",color = colorRampPalette(paletteer::paletteer_d("RColorBrewer::PuBu"))(100),#c("white", "#FF6666") paletteer::paletteer_d("RColorBrewer::Purples")
breaks=(c(seq(0, 0.1, length.out = 5),seq(0.11, 0.2, length.out = 5),seq(0.21, 0.3, length.out = 20),
seq(0.31, 0.4, length.out = 30),seq(0.41, 1, length.out = 40))),
file="Fig1E_Bray-Curtis distance.pdf")



6.Fig S1F(correlations)
####Top expression gene
#correlation among all single-cell RNA expression
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
#celltype2
Idents(combined)<- "celltype2"
av.exp<- AverageExpression(combined)$RNA
# av.exp<- av.exp[which(row.names(av.exp)%in% features),]

features=names(tail(sort(apply(av.exp, 1, sd)),2000))
av.exp<- av.exp[which(row.names(av.exp)%in% features),]
cor.av.exp <- cor(av.exp, method= "spearman")
##change cell order
cor.av.exp <- cor.av.exp[,c(16,1,14,11,13,7,5,12,8,9,10,6,4,15,2,3)] 
cor.av.exp <- cor.av.exp[c(16,1,14,11,13,7,5,12,8,9,10,6,4,15,2,3),]

pheatmap(cor.av.exp,color = colorRampPalette(c("white","red3"))(10),
breaks=(c(seq(0.41, 0.7, length.out = 4),seq(0.71, 1, length.out = 6))),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = TRUE,       
border_color = NA,   scale = "none",      
height = 6,width=7,
file="spearman_cor_RNA.pdf")


##ATAC spearman corrlation
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
##call peak by celltype2
DefaultAssay(combined) <- 'ATAC'
# call peaks for each celltype in myeloid using MACS2
peaks <- CallPeaks(combined, group.by = "celltype2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks1 <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks1 <- subsetByOverlaps(x = peaks1, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(combined),
  features = peaks1,
  cells = colnames(combined)
)

fragpath <- Fragments(combined)
# create a new assay using the MACS2 peak set and add it to the Seurat object
combined[["Peaks1"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(combined,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
#cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/
##cor
Idents(combined) <- "celltype2"
cell=combined$celltype2
cell_u=c("pDC","CM","ncM","cDC","B plasma","B naive","B memory",
         "NK","CD8 Tem","CD8 Tcm","CD8 MAIT",
         "CD8 Treg","CD8 naive","CD4 Treg","CD4 naive","CD4 Tm")
DefaultAssay(combined) <- 'Peaks1'
peakMat=GetAssayData(combined,slot="data")

PEAK=c()
for(i in 1:length(cell_u)){
    print(i)
  peakMat1=peakMat[,which(cell==cell_u[i])]
  peakMat2=apply(peakMat1,1,mean)
  PEAK=rbind(PEAK,peakMat2)
}
PEAK=t(PEAK)
colnames(PEAK)=cell_u
R_PEAK=cor(PEAK,method="spearman")
write.table(R_PEAK,"spearman_cor_peak.txt",quote=F,sep="\t")
pheatmap(R_PEAK,color = colorRampPalette(c("white","Orange2"))(10),
#breaks=(c(seq(0.51, 0.8, length.out = 2),seq(0.81, 1, length.out = 8))),
cluster_rows = FALSE,   cluster_cols = FALSE,      
show_rownames = TRUE,show_colnames = TRUE,       
border_color = NA,   scale = "none",      
height = 6,width=7,
file="spearman_cor_PEAK.pdf")





