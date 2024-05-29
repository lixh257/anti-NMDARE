
Part3:
1. Reduction and clustering
2.Annotation


###
1. Reduction and clustering
library(Seurat)
library(ggplot2)
library(Signac)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
combined=readRDS("Merge_after_filter1.rds")
combined.list <- SplitObject(combined, split.by = "dataset")
for (i in 1:length(combined.list )) {
 DefaultAssay(combined.list [[i]]) <- "RNA"
combined.list [[i]] <- SCTransform(combined.list [[i]], verbose = FALSE)
}
for (i in seq_len(length(combined.list ))) {
  DefaultAssay(combined.list [[i]]) <- "SCT"
}
combined.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features)
#integrate RNA using rpca
combined_list <- lapply(
  X = combined.list,
  FUN = ScaleData,
  features = combined.features,
  verbose = FALSE
)
combined_list <- lapply(
  X =combined.list,
  FUN = RunPCA,
  features = combined.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list =combined_list,
  normalization.method = "SCT",
  anchor.features = combined.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)
combined1<- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA",
  dims = 1:30
)
#combined <- RunHarmony(combined, group.by.vars ="dataset",reduction.save="harmony",assay.use ="SCT") 
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(combined1) <- "ATAC"
combined1<- RunTFIDF(combined1)
combined1<- FindTopFeatures(combined1, min.cutoff = "q25")
combined1 <- RunSVD(combined1)
#integrate embeddings and output new object to prevent overwriting integrated RNA
combined1_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI",
  reductions = combined1@reductions$lsi
)
#copy integrated LSI from duplicate seurat object to original object
combined1@reductions$integratedLSI <- combined1_atac@reductions$integratedLSI
# assess the correlation between each LSI component and sequencing depth using the DepthCor() function
pdf("DepthCor.pdf")
DepthCor(combined1,n=50,reduction = 'integratedLSI')
dev.off()

DefaultAssay(combined1) <- "ATAC" 
combined1 <- RunUMAP(combined1, reduction = 'integratedLSI', dims = 2:50, reduction.name = "inter.umap.atac", reduction.key = "interatacUMAP_")
pdf("atac_umap_test.pdf",width=10, height=10)
DimPlot(combined1, reduction = "inter.umap.atac", group.by = "Label", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC") 
dev.off()

DefaultAssay(combined1) <- "integratedRNA" 
combined1 <- RunPCA(combined1) #reductions names "pca"
#ElbowPlot(combined1 ,  ndims = 50)
# Printing out the most variable genes driving PCs 
pdf("pc_heatmap.pdf",width=9,height=20)
DimHeatmap(combined1,reduction = 'pca', dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
combined1 <- RunUMAP(combined1,reduction = 'pca',dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
pdf("rna_umap_test.pdf",width=10, height=10)
DimPlot(combined1, reduction = "umap.rna", group.by = "Label", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA") 
dev.off()

combined1<- FindMultiModalNeighbors(combined1, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:30, 2:50))
combined1 <- RunUMAP(combined1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pdf("WNN_umap_test.pdf",width=10, height=10)
DimPlot(combined1, reduction = "wnn.umap", group.by = "Label", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
dev.off()


combined1<- FindClusters(combined1, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(combined1$seurat_clusters)
pdf("WNN_umap_test1.pdf",width=10, height=10)
DimPlot(combined1, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
dev.off()


setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/")
p1 <- DimPlot(combined1, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA") 
p2 <- DimPlot(combined1, reduction = "inter.umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC") 
p3 <- DimPlot(combined1, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
pdf("Merge_doimplot_all_clusters.pdf",width=30, height=10)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()


p1 <- DimPlot(combined1, reduction = "umap.rna", group.by = "dataset", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA") 
p2 <- DimPlot(combined1, reduction = "inter.umap.atac", group.by = "dataset", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC") 
p3 <- DimPlot(combined1, reduction = "wnn.umap", group.by = "dataset", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
pdf("Merge_doimplot_all_dataset.pdf",width=30, height=10)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()
saveRDS(combined1,"Merge_after_filter1.rds")


##other plot
combined=combined1
pdf("Merge_doimplot_all_AE_dataset.pdf",width=15,height=5)
DimPlot(combined1, reduction = "wnn.umap", group.by = 'AE', pt.size = 0.1,split.by = 'AE')
dev.off()

pdf("Merge_doimplot_all_AE.pdf",width=10,height=5)
DimPlot(combined, reduction = "wnn.umap", group.by = 'AE', pt.size = 0.1,split.by = 'AE')
dev.off()

pdf("Merge_doimplot_all_Sex_dataset.pdf",width=15,height=5)
DimPlot(combined, reduction = "wnn.umap", group.by = 'Sex', pt.size = 0.1,split.by = 'dataset')
dev.off()
pdf("Merge_doimplot_all_Sex.pdf",width=10,height=5)
DimPlot(combined, reduction = "wnn.umap", group.by = 'Sex', pt.size = 0.1,split.by = 'Sex')
dev.off()

pdf("Merge_doimplot_all_Label_dataset.pdf",width=15,height=5)
DimPlot(combined, reduction = "wnn.umap", group.by = 'Label', pt.size = 0.1,split.by = 'dataset')
dev.off()
pdf("Merge_doimplot_all_Label.pdf",width=6,height=5)
DimPlot(combined, reduction = "wnn.umap", group.by = 'Label', pt.size = 0.1)
dev.off()



2.Annotation
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/")
library(Seurat)
library(dplyr)
combined=readRDS("Merge_after_filter1.rds")
pdf("Merge_doimplot_all_clusters1.pdf",width=12, height=10)
DimPlot(combined, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
dev.off()
pdf("Merge_doimplot_all_test1.pdf",width=12, height=10)
DimPlot(combined, reduction = "wnn.umap", group.by = "celltype2", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN") 
dev.off()

Idents(combined)="seurat_clusters"
DefaultAssay(combined) <- "RNA" 
pbmc.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
pbmc.markers %>%     
group_by(cluster) %>%     
top_n(n = 10, wt = avg_log2FC) -> top10 
write.table(top10,"markers_top10.txt",quote=F,sep="\t",row.names=F)


##大类：T cell(CD4,CD8),NK,B,Myeloid，其他可能没筛干净的血细胞等（血小板platelet，浆细胞plasma cell，红细胞Erythroid，祖细胞Early Progenitor）

markerGenes  <- c(
##1.other cell
 #1.1 Progenitor/stem cell祖细胞干细胞:
"CD34","PRSS57","SOX4",
#1.2 Early erythrocyte/erythrocyte-like,有核的红细胞:
"HBD","GATA1","KIT","TFRC","PIEZO1","CGA",
##2.immune cell
##2.1 Bcell:
"MS4A1","CD79A","BLNK","CD19","CD79B",
#2.1.1 B naive:
"CD27","CD38","IgD","CD72","SELL","TCL1A",
#2.1.2 B memory:
"CD27","BANK1","CD79A","MS4A1","ITGB2",
#2.1.3 B plasma
"CD38","CD19","MZB1",
##2.2 NKT cell
"CD3D","CD3E",
#2.2.1 NK
"NKG7","GNLY","KLRB1","KLRF1","KLRD1","GZMK","GZMB","PRF1",
#2.2.2 CD4+Tn
"CD4","CCR7","SELL","LEF1",
#2.2.3 CD8+Tn
"CD8A","CD8B","CCR7","SELL","LEF1",
#2.2.4 CD4+Tm
"CD4","IL7R","TCF7",
#2.2.5 CD8+Tcm(central memory)
"CD3D","CD8B","NKG7","GZMK","GZMH",
#2.2.6 C8+Tem(effector memory)
"GZMB","GZMH","PRF1","NKG7",
#2.2.7 Treg
"FOXP3","RTKN2",
#2.2.8 gdT(gama delta T cell)
"TRDV2","TRGV9",
#2.2.9 MAIT
"SLC4A10","IL7R","KLRB1",
##2.3 Myeloid
"LYZ","CD14","HLA-DQA1","S100A9","VCAN","FCGR3A",
#2.3.1 Megakaryocytes巨核细胞
"GNG11","ITGA2B",
#2.3.2 Classical Monocytes
"VCAN","S100A9","CD14","ITGAM","ITGAX","IL3RA",
#2.3.3 Non-classical Monocytes
"CD68","S100A12","MTSS1","SERPINA1","MS4A7","FCGR3A",
#2.3.4 cDC 
"FLT3","FCER1A","CLEC10A","DAPP1",
#2.3.5 pDC
"IRF8","IRF7","LILRA4",
#2.3.6 Macrophage(可能没有)
"CD163","CCR7","IL7R","CD14","CD16"
)
markerGenes=unique(markerGenes)
#IgD, CD16(FCGR3A), CD41a(ITGA2B), CD11b(ITGAM), CD11c(ITGAX), CD123(IL3RA), CD127(IL7R)


mainmarker=c(
"TCL1A","BANK1","PAX5","MS4A1","CD79A","MZB1","IGHM",#B-Cell
"CD3D","CD3E","CD3G","CD8A", "CD4", "IL7R","CCL5","CCR6","CCR7","RTKN2",#TCells
 "NKG7","GNLY","KLRB1","KLRF1","KLRD1",#NK
"LYZ","CD14","CEBPB","IRF7","LILRA4", "FCER1A","CLEC10A","HLA-DQA1","DAPP1","S100A9","VCAN","MTSS1","SERPINA1","MS4A7","FCGR3A"#Myeloid
)
Idents(combined)="seurat_clusters"
pdf("Annotation_Dotplot_main.pdf",height=10,width=15)
DotPlot(combined, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()


8,9,16-B
0,12,15,19-Myeloid
1:7,10:11,13,14,17,18,20-TNK

Idents(combined)=combined$seurat_clusters
combined <- RenameIdents(combined, '8' = 'B','9' = 'B','16' = 'B') 
combined <- RenameIdents(combined, '0' = 'Myeloid', '12' ='Myeloid', '15' = 'Myeloid', '19' = 'Myeloid') 
combined <- RenameIdents(combined, '1' = 'TNK', '2' = 'TNK', '3' = 'TNK', '4' = 'TNK', '5' = 'TNK', '6' = 'TNK',
                                   '7' = 'TNK', '10' = 'TNK', '11' = 'TNK', '13' = 'TNK', '14' = 'TNK','17' = 'TNK',
                                   '18' = 'TNK','20' = 'TNK') 
combined$celltype1=Idents(combined)
pdf("Merge_doimplot_all_Celltype1.pdf",width=10,height=8)
DimPlot(combined, reduction = "wnn.umap", group.by = "celltype1", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("celltype1") 
dev.off()
saveRDS(combined,"Merge_after_filter1.rds")

setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/")
combined=readRDS("Merge_after_filter1.rds")

2.2 Sub Cell
2.2.1 B Cell
8,9,16-B
0,12,15,19-Myelo
1:7,10:11,13,14,17,18,20-TNK

combined_B=subset(x=combined,
celltype1 == "B"
) 
DefaultAssay(combined_B) <- "RNA"

combined_B <- FindMultiModalNeighbors(combined_B, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:40, 2:50))
combined_B <- RunUMAP(combined_B, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combined_B <- FindClusters(combined_B, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(combined_B$seurat_clusters)

markerGenes  <- c(
"CD19", #B
"CD27","CD38","SDC1","MZB1",#Bplasma
"CD27","CD38","CD72","SELL","TCL1A",#Bnaive
"CD27","BANK1","CD79A","MS4A1","ITGB2"#Bmemory
)

Idents(combined_B)="seurat_clusters"
markerGenes=unique(markerGenes)
pdf("Dotplot_B_clusters.pdf",width=6,height=4)
DotPlot(combined_B, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()

pdf("Dimplot_B_clusters.pdf",height=5,width=6)
DimPlot(combined_B, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()

###5-Bp;0,3,7,9-Bn,1,2,4,6,8-Bm

Idents(combined_B)=combined_B$seurat_clusters
combined_B <- RenameIdents(combined_B, '5' = 'B plasma')
combined_B <- RenameIdents(combined_B, '0' = 'B naive', '3' = 'B naive', '7' = 'B naive', '9' = 'B naive') 
combined_B <- RenameIdents(combined_B, '2' = 'B memory','1' = 'B memory','4' = 'B memory','6' = 'B memory','8' = 'B memory') 
combined_B$celltype2=Idents(combined_B)
pdf("Dimplot_B_celltype2.pdf",height=5,width=6)
DimPlot(combined_B, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()

2.2.2 T Cell
combined_T=subset(x=combined,
celltype1 == "TNK"
) 
DefaultAssay(combined_T) <- "RNA"
combined_T<- FindMultiModalNeighbors(combined_T, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:40, 2:50))
combined_T<- FindClusters(combined_T, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.8)
combined_T <- RunUMAP(combined_T, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
table(combined_T$seurat_clusters)

pdf("Dimplot_T_clusters.pdf",width=5,height=5)
DimPlot(combined_T, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()


markerGenes=c(
##2.2 NKT cell
"CD3D","CD3E",
#2.2.1 NK
"NKG7","GNLY","KLRB1","KLRF1","KLRD1","GZMK","GZMB","PRF1",
#2.2.2 CD4+Tn
"CD4","CCR7","SELL","LEF1",
#2.2.3 CD8+Tn
"CD8A","CD8B","CCR7","SELL","LEF1",
#2.2.4 CD4+Tm
"CD4","IL7R","TCF7",
#2.2.5 CD8+Tcm(central memory)
"CD3D","CD8B","NKG7","GZMK","GZMH",
#2.2.6 C8+Tem(effector memory)
"GZMB","GZMH","PRF1","NKG7",
#2.2.7 Treg
"FOXP3","RTKN2",
#2.2.8 gdT(gama delta T cell)
"TRDV2","TRGV9",
#2.2.9 MAIT
"SLC4A10","IL7R","KLRB1"
)
markerGenes=unique(markerGenes)
pdf("Dotplot_T_clusters.pdf",width=8,height=6)
DotPlot(combined_T, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()

pdf("Dotplot_T_clusters_1.pdf",width=8,height=6)
DotPlot(combined_T, features = markerGenes[c(8:10,15,16,19,22:24)], cols = c("white", "blue"), idents=c(2,3,7,8,15,19,20,23,24,25),dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()
pdf("Dotplot_T_clusters_2.pdf",width=8,height=6)
DotPlot(combined_T, features = markerGenes, cols = c("white", "blue"), idents=c(2,6,8,7,9,15,19,21,24,25),dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()


12,17-MAIT;   4,22-NK; 
14,1,10,11-CD4naive; 6,5,9,-CD4Tm, 18,16-Treg;
21-CD8Treg;13,0-CD8naive;8,15,20,23,25-CD8Tcm;2,3,7,19,24-CD8Tem;

Idents(combined_T)=combined_T$seurat_clusters
combined_T <- RenameIdents(combined_T, '22' = 'NK', '4' = 'NK',  
                                       '16'='CD4 Treg','18'='CD4 Treg','21'='CD8 Treg',  
                                       '12'="CD8 MAIT",'17'="CD8 MAIT")
combined_T <- RenameIdents(combined_T, '14' = 'CD4 naive','1' = 'CD4 naive', '10' = 'CD4 naive','11' = 'CD4 naive',
                                       '13' = 'CD8 naive','0' = 'CD8 naive') 
combined_T <- RenameIdents(combined_T, '6' = 'CD4 Tm','5' = 'CD4 Tm', '9' = 'CD4 Tm') 
combined_T <- RenameIdents(combined_T, '2' = 'CD8 Tem','3' = 'CD8 Tem','7' = 'CD8 Tem', 
                                       '19' = 'CD8 Tem','24' = 'CD8 Tem', 
                                       '8' = 'CD8 Tcm','15' = 'CD8 Tcm', '20' = 'CD8 Tcm', '23' = 'CD8 Tcm', '25' = 'CD8 Tcm') 

combined_T$celltype2=Idents(combined_T)
pdf("Dimplot_T_celltype.pdf",height=5,width=6)
DimPlot(combined_T, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()
saveRDS(combined_T,"Merge_5lib_TNK.rds")
saveRDS(combined_B,"Merge_5lib_Bcell.rds")


2.2.3 Myeloid
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb) 
library(EnsDb.Hsapiens.v86) 
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(scDblFinder)
library(dittoSeq)
combined_M=subset(x=combined,celltype1=="Myeloid")
DefaultAssay(combined_M) <- "RNA"
combined_M<- FindMultiModalNeighbors(combined_M, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:30, 2:50))
combined_M <- RunUMAP(combined_M, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combined_M<- FindClusters(combined_M, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(combined_M$seurat_clusters)

pdf("Dimplot_M_clusters.pdf",width=5,height=5)
DimPlot(combined_M, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black")
dev.off()
pdf("Dimplot_M_clusters1.pdf",width=5,height=5)
DimPlot(combined_M, reduction = "wnn.umap", group.by = 'Label', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black")
dev.off()

pdf("Dimplot_M_clusters_1.pdf",width=25,height=5)
DimPlot(combined_M, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black",split.by="Label")
dev.off()
DefaultAssay(combined_M)="RNA"
markerGenes  <- c(
##2.3 Myeloid
#"LYZ","CD14","HLA-DQA1","S100A9","VCAN","FCGR3A",
#2.3.1 Megakaryocytes巨核细胞
"GNG11","ITGA2B",
#2.3.2 Classical Monocytes
"CD14","ITGAM","ITGAX","IL3RA",
#2.3.3 Non-classical Monocytes
"CD68","S100A12","MTSS1","SERPINA1","MS4A7","FCGR3A",
#2.3.4 cDC 
"FLT3","FCER1A","CLEC10A","DAPP1",
"CLEC9A",
#2.3.5 pDC
"IRF8","IRF7","LILRA4",
#2.3.6 Macrophage(可能没有)
"CD163","CCR7","IL7R","CD14","FCGR3A"
)
markerGenes  <- c(
"CD14","VCAN",
"MS4A7","FCGR3A",
"FLT3","FCER1A","CLEC10A","CLEC9A",
"IRF8","IRF7","LILRA4","CLEC4C",
"GNG11","ITGA2B","PF4",
"CSF1R","CD209","MRC1",
##Macro
"CD163","CCR2","CCL2",#(PMID:32561858
"INHBA","IL1RN","CCL4","NLRP3","EREG","IL1B",
"LYVE1","PLTP","SEPP1","C1QC","C1QA","APOE"#Macro(,PMID:33545035))
#"C1QB","APOE"
)
#MRC1(CD206),
Idents(combined_M)="seurat_clusters"
markerGenes=unique(markerGenes)
pdf("Dotplot_M_clusters.pdf",width=8,height=6)
DotPlot(combined_M, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()



pbmc.markers <- FindAllMarkers(combined_M, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
pbmc.markers %>%     
group_by(cluster) %>%     
top_n(n = 10, wt = avg_log2FC) -> top10 
write.table(pbmc.markers ,"combined_M_markers_markers .txt",quote=F,sep="\t",row.names=F)


markerGenes=c(
"IRF4","IRF7","IRF8","SPIB","SOX4","LILRA4","SLC15A4","PLD4","CCDC50","IL3RA","LY9", "SELL","GAS6",#pDC
"HMGA1","PFDN1","IRF4","CD1C","CLEC10A","FCGR2B","ADAM8","FCER1A","AXL","ADAM28","LY86","TIMM13","ARAF",#cDC2
"BATF3","ID2","ETV3","XCR1","CLEC9A","PTDSS1","SCARB1","IL6ST","CD40","TNFRSF10B","IDO1","CST7","CLIC2","NET1","ANXA6",#cDC1
"CEBPD","FCN1","CD14","CD36","SELL","S100A8","S100A12","CLEC12A","MS4A6A","CXCL14",##class Mono
"TCF7L2","KLF3","IKZF1","FLI1","FCN1","FCGR3A","FCGR3B","CX3CR1","IFITM1","ICAM2","MTSS1","CDKN1C","CDH23","SLC44A2",#non-class mono
"CEBPD","ZFP36L2","FCN1","CD14","CD300E",#inter#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10218702/
"TNF","FCGR1A","IRF1","HLA-DPB1","CD86","MARCO","IL2RA",##M1 macro
"CXCR4","IL27RA","CSF1R","CTSD","HMOX1","PPARG","LIPA","CLEC7A","F13A1","MAF","MS4A4A" #M2 macro
)

markerGenes=unique(markerGenes)
pdf("Dotplot_M_clusters_1.pdf",width=15,height=6)
DotPlot(combined_M, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, cluster.idents  =TRUE) + RotatedAxis()
dev.off()

##CM,ncM,inter_M:#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10218702/
markerGenes=c(
"CD14","FCGR3A",
"FCGR1A","CD86","TNFRSF1A","TNFRSF1B","CD38",#CD64-FCGR1A,TNFR1-TNFRSF1A,TNFR2-TNFRSF1B,HLADR-CD38
"CD163","CD36",
"CCR2","CCR5","CX3CR1","ITGAM","SELL",#CD11b-ITGAM,CD62L-SELL
"IL6","CCL2","CCL3","IL1B"
)
pdf("Dotplot_M_clusters_2.pdf",width=8,height=4)
DotPlot(combined_M, features = markerGenes, cols = c("white", "blue"), dot.scale = 8, 
        cluster.idents  =TRUE,ident=c(0:7)) + RotatedAxis()
dev.off()



DefaultAssay(combined_M)="RNA"
library(dittoSeq)
pdf("Featureplot_M_clusters_2.pdf",width=12,height=6)
dittoDimPlot(combined_M, 
             var = c("TNFRSF1A","CCR5","CX3CR1"),#,"CD14","FCGR3A","CD86"
             size = 1,
             theme = theme_classic(),
             min.color = "white smoke",
             max.color = "blue",
             main = "Markers of Myeloid",
           #  min = 0,
            # max = 5,
             order = "increasing",
            reduction = "wnn.umap")#这个函数的优势
dev.off()


9-pDC,8-cDC,6-Non-classical Monocytes ;
0:5,7-Classical Monocytes
Idents(combined_M)=combined_M$seurat_clusters
combined_M<- RenameIdents(combined_M, '6' = 'ncM',  '8' = 'cDC' ,'9'='pDC')
combined_M <- RenameIdents(combined_M,'0'="CM",'1'="CM",'2'="CM",
                                      '4'="CM",'3'="CM",'5'="CM",'7'="CM")
combined_M$celltype2=Idents(combined_M)
pdf("Dimplot_M_celltype2.pdf",height=5,width=8)
DimPlot(combined_M, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()

saveRDS(combined_M,"Merge_5lib_Myeloid.rds")


2.3 Rename
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/")
combined_M=readRDS("Merge_5lib_Myeloid.rds")
combined_T=readRDS("Merge_5lib_TNK.rds")
combined_B=readRDS("Merge_5lib_Bcell.rds")
combined=readRDS("Merge_after_filter1.rds")
unique(combined_M$celltype2)
unique(combined_T$celltype2)
unique(combined_B$celltype2)
[1] CM  cDC ncM pDC
[1] CD4 naive CD4 Tm    CD8 naive CD8 Treg  CD8 Tem   CD8 Tcm   CD8 MAIT NK        CD4 Treg 
[1] B memory B naive  B plasma
combined$celltype2=0
combined@meta.data[rownames(combined_M@meta.data)[which(combined_M$celltype2=='ncM')],60]='ncM'
combined@meta.data[rownames(combined_M@meta.data)[which(combined_M$celltype2=='CM')],60]='CM'
combined@meta.data[rownames(combined_M@meta.data)[which(combined_M$celltype2=='cDC')],60]='cDC'
combined@meta.data[rownames(combined_M@meta.data)[which(combined_M$celltype2=='pDC')],60]='pDC'

combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD4 Tm')],60]='CD4 Tm'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD8 Tem')],60]='CD8 Tem'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD8 Tcm')],60]='CD8 Tcm'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD4 naive')],60]='CD4 naive'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD8 naive')],60]='CD8 naive'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='NK')],60]='NK'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD4 Treg')],60]='CD4 Treg'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD8 Treg')],60]='CD8 Treg'
combined@meta.data[rownames(combined_T@meta.data)[which(combined_T$celltype2=='CD8 MAIT')],60]='CD8 MAIT'
#B memory  B naive B plasma
combined@meta.data[rownames(combined_B@meta.data)[which(combined_B$celltype2=='B memory')],60]='B memory'
combined@meta.data[rownames(combined_B@meta.data)[which(combined_B$celltype2=='B naive')],60]='B naive'
combined@meta.data[rownames(combined_B@meta.data)[which(combined_B$celltype2=='B plasma')],60]='B plasma'
table(combined$celltype2)

pdf("Dimplot_all_celltype.pdf",height=5,width=6)
DimPlot(combined, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()

#Save
saveRDS(combined,"Merge_after_filter1.rds")

#cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_after_filter1.rds /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/

