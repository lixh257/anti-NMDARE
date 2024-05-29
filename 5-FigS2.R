
Fig S2
#####################
library(Seurat)
library(ggplot2)
library(dittoSeq)

setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
TCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
Myeloid=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")

marker_gene_B=c("TCL1A","MS4A1","MZB1","SDC1","CD27","CD38")
marker_gene_T=c("CD3D","CD8A","IL7R","CCR7","GZMB","GZMK","KLRB1","GNLY","FOXP3")
marker_gene_M=c("CD14","FCGR3A","FCER1A","GAS6")

sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
            "#2DB7A3","#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

1.Fig S2A
DefaultAssay(BCell)="SCT"#"SCT"/intergrate_RNA
pdf("FigS2_B_marker_Featureplot_RNA_TCL1A.pdf",width=4,height=3)
dittoDimPlot(BCell, 
             marker_gene_B[1],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_B_marker_Featureplot_RNA_MS4A1.pdf",width=4,height=3)
dittoDimPlot(BCell, 
             marker_gene_B[2],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_B_marker_Featureplot_RNA_MZB1.pdf",width=4,height=3)
dittoDimPlot(BCell, 
             marker_gene_B[3],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_B_marker_Featureplot_RNA_SDC1.pdf",width=4,height=3)
dittoDimPlot(BCell, 
             marker_gene_B[4],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()


2.Fig S2B
marker_gene_T=c("CD4","CD8A","IL7R","CCR7","GZMB","GZMK","KLRB1","GNLY","FOXP3")
DefaultAssay(TCell)="SCT"#"SCT"/intergrate_RNA
pdf("FigS2_T_marker_Featureplot_RNA_CD4.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[1],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_T_marker_Featureplot_RNA_CD8A.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[2],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_T_marker_Featureplot_RNA_CCR7.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[4],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
            #max=8,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_T_marker_Featureplot_RNA_IL7R.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[3],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()


marker_gene_T=c("CD4","CD8A","CCR7","IL7R","FOXP3","RTKN2","GZMB","GZMK","KLRB1","GNLY")
pdf("FigS2_T_marker_Featureplot_RNA_FOXP3.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[5],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_T_marker_Featureplot_RNA_RTKN2.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[6],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

marker_gene_T=c("CD4","CD8A","CCR7","IL7R","FOXP3","RTKN2",
                "GZMB","GZMK","GZMH","KLRB1","NCR3","GNLY")
pdf("FigS2_T_marker_Featureplot_RNA_GZMB.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[7],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_T_marker_Featureplot_RNA_GZMK.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[8],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_T_marker_Featureplot_RNA_GZMH.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[9],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_T_marker_Featureplot_RNA_KLRB1.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[10],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_T_marker_Featureplot_RNA_NCR3.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[11],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_T_marker_Featureplot_RNA_GNLY.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             marker_gene_T[12],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_T_marker_Featureplot_RNA_NKG7.pdf",width=4,height=3)
dittoDimPlot(TCell, 
             "NKG7",
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

3.Fig S2C
Myeloid=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
marker_gene_M=c("CD14","FCGR3A","FCER1A","GAS6")
DefaultAssay(Myeloid)="SCT"
pdf("FigS2_M_marker_Featureplot_RNA_CD14.pdf",width=4,height=3)
dittoDimPlot(Myeloid, 
             marker_gene_M[1],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_M_marker_Featureplot_RNA_FCGR3A.pdf",width=4,height=3)
dittoDimPlot(Myeloid, 
             marker_gene_M[2],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
pdf("FigS2_M_marker_Featureplot_RNA_FCER1A.pdf",width=4,height=3)
dittoDimPlot(Myeloid, 
             marker_gene_M[3],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()

pdf("FigS2_M_marker_Featureplot_RNA_GAS6.pdf",width=4,height=3)
dittoDimPlot(Myeloid, 
             marker_gene_M[4],
             size = 1,
             #ncol=1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             #max=1.5,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()