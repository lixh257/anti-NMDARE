##Fig1:Annotation result and cell distribution

1.Fig1A. Workflow by AI

2. Fig 1B (UMAP of celltypes)
library(Seurat)
library(ggplot2)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_after_filter1.rds")

pal <- paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20, 3,13,19, 7,15,17,5,14,11, 12,16,18,8)]
Idents(combined)="celltype2"
cell=c("B naive","B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg","CD8 naive",
"CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK","CM","ncM","cDC","pDC")
names(cell) <- levels(combined)
combined <- RenameIdents(combined, cell)
levels(combined) <- c("B naive","B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg","CD8 naive",
"CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK","CM","ncM","cDC","pDC")
p1 <- DimPlot(combined, reduction = "umap.rna", cols= pal, 
group.by = "celltype2", label = TRUE, label.size = 5, 
pt.size = 0.5,repel = TRUE) + ggtitle("RNA")+ 
NoLegend() +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
p2 <- DimPlot(combined, reduction = "inter.umap.atac", cols= pal, 
group.by = "celltype2", label = TRUE, label.size = 5, repel = TRUE,pt.size = 0.5) + 
ggtitle("ATAC")+
 NoLegend() +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
p3 <- DimPlot(combined, reduction = "wnn.umap", cols= pal, 
group.by = "celltype2", label = TRUE, label.size = 5, 
repel = TRUE,pt.size = 0.5) + 
ggtitle("WNN") +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框

pdf("Fig1B_Merge_doimplot_all_celltype2.pdf",width=35, height=10)
p1 + p2 + p3  & theme(plot.title = element_text(hjust = 0.5))
dev.off()



3.Fig 1C (marker gene expression)
markerGenes  <- c(
##1.other cell
#"CD34","SOX4", #1.1 Progenitor/stem cell祖细胞干细胞:
##2.immune cell
##2.1 Bcell:
"TCL1A","CD79A","MS4A1","MZB1",
"CD27","CD38","SDC1",#Bplasma
##2.2 NKT cell
"CD3E","CD3D",
"CD4","CCR7",  #2.2.2 CD4+Tn
"CD4","IL7R",#2.2.4 CD4+Tm
"FOXP3","RTKN2",#2.2.7 Treg
"CD8A",#2.2.3 CD8+Tn
"GZMK","GZMH",#2.2.5 CD8+Tcm(central memory)
"GZMB","GZMH",#2.2.6 C8+Tem(effector memory)
"IL7R","KLRB1","NCR3",#2.2.9 MAIT
"NKG7","GNLY",#2.2.1 NK
##2.3 Myeloid
#"ITGA2B","PF4",#2.3.1 Megakaryocytes巨核细胞
"VCAN","LYZ","CD14",#2.3.2 Classical Monocytes
"MS4A7","FCGR3A",#2.3.3 Non-classical Monocytes
"FLT3","FCER1A","CLEC10A",#2.3.4 cDC 
"GAS6","IRF8","IRF7","LILRA4"#2.3.5 pDC
)
markerGenes=unique(markerGenes)
#dotplot的改造  
DefaultAssay(combined) <- "RNA"
Idents(combined)="celltype2"
pdf("Fig1C_Dotplot_Marker11_1.pdf",width=7,height=8)
DotPlot(combined, features = markerGenes,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,10, 20, 30, 50,75),range = c(-1,7))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)
dev.off()



4.Fig 1D (marker gene’s Chromatin accessibility)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_after_filter1.rds")
Idents(combined)="celltype2"
DefaultAssay(combined) <- 'Peaks1'##call peak by celltype in FigS1.R

gene.activities <- GeneActivity(combined)
# add gene activities as a new assay
combined[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
#combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/Merge_after_filter.rds")
###get motif & add chromvar
library(chromVARmotifs)
library(dplyr)
library(tidyr)
library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
)

motif.mat <- CreateMotifMatrix(
  features = StringToGRanges(rownames(combined)),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

combined <- RegionStats(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

#RunChromVAR to caculate the motif activities
combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  new.assay.name = "chromvar"
)
saveRDS(combined,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_after_filter1.rds")


library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_after_filter1.rds")
markerGenes  <- c(
##1.other cell
#"CD34","SOX4", #1.1 Progenitor/stem cell祖细胞干细胞:
##2.immune cell
##2.1 Bcell:
"TCL1A","CD79A","MS4A1","MZB1",
"CD27","CD38","SDC1",#Bplasma
##2.2 NKT cell
"CD3E","CD3D",
"CD4","CCR7",  #2.2.2 CD4+Tn
"CD4","IL7R",#2.2.4 CD4+Tm
"FOXP3","RTKN2",#2.2.7 Treg
"CD8A",#2.2.3 CD8+Tn
"GZMK","GZMH",#2.2.5 CD8+Tcm(central memory)
"GZMB","GZMH",#2.2.6 C8+Tem(effector memory)
"IL7R","KLRB1","NCR3",#2.2.9 MAIT
"NKG7","GNLY",#2.2.1 NK
##2.3 Myeloid
#"ITGA2B","PF4",#2.3.1 Megakaryocytes巨核细胞
"VCAN","LYZ","CD14",#2.3.2 Classical Monocytes
"MS4A7","FCGR3A",#2.3.3 Non-classical Monocytes
"FLT3","FCER1A","CLEC10A",#2.3.4 cDC 
"GAS6","IRF8","IRF7","LILRA4"#2.3.5 pDC
)
markerGenes <- unique(markerGenes)

########Marker gene Track plot
DefaultAssay(combined)="Peaks1"
cell=c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
names(cell) <- levels(combined)
combined <- RenameIdents(combined, cell)
combined <- RegionStats(combined, genome = BSgenome.Hsapiens.UCSC.hg38)
Gene=markerGenes

combined <- LinkPeaks(
  object = combined,
  peak.assay = "Peaks1",
  expression.assay = "SCT",
  genes.use = Gene,
  pvalue_cutoff = 0.01,
  score_cutoff = 0
  )
Idents(combined) <- "celltype2"
DefaultAssay(combined)="Peaks1"

levels(combined) <- rev(c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC"))
#pdf("Fig1_MarkerGene_track_1000pb_1.pdf",height=6,width=8)
#for(i in c(1:34)){
  i=1
pdf(paste0("Fig1_MarkerGene_track_",Gene[i],".pdf"),height=6,width=3)
p=CoveragePlot(
   object = combined,
   region = Gene[i], 
   #features = Gene[i],
   extend.upstream = -33000,###change this value by every gene's feature
   extend.downstream =-7000,###change this value by every gene's feature
   #ymax = 120,##change normalized signal range
   peaks = FALSE,
   links=FALSE,
   idents=factor(cell,levels=cell),
   window=300
   #col=c('#E5D2DD','#D6E7A3'),
   #scale.factor = 1000000##change signal showing scale
   #heights = c(5, 5,  2, 5),
   #widths = c(12, 2)
 )+ scale_fill_brewer(type = "qual", palette = 1)
 print(p)
 dev.off()
#}
#dev.off()


5.Fig 1E  (UAMP the cell with different Samples)
library(Seurat)
library(ggplot2)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter2.rds")
pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#539092","#8BCFCC","#51ADCF","#0278AE","#01A9B4","#A3D8F4")
pal_other=c("#D8D9DA")
combined$Patient1=rep("Other",length(combined$Label))
combined$Patient1[which(combined$Label=="Patient1")]="Patient1"
Idents(combined)="Patient1"
pdf("Fig1E_Pat1.pdf",width=6,height=5)
p=DimPlot(combined, reduction = "wnn.umap", 
cols=c(alpha(pal_other[1],0.5),alpha(pal_pat[1],0.8)),
group.by = "Patient1", label = TRUE, label.size = 1, 
order=c("Patient1","Other"),
repel = TRUE,pt.size = 0.05) + 
ggtitle("WNN") +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
print(p)
dev.off()


combined$Patient3=rep("Other",length(combined$Label))
combined$Patient3[which(combined$Label=="Patient3")]="Patient3"
Idents(combined)="Patient3"
pdf("Fig1E_Pat3.pdf",width=6,height=5)
p=DimPlot(combined, reduction = "wnn.umap", 
cols=c(alpha(pal_other[1],0.5),alpha(pal_pat[2],0.8)),
group.by = "Patient3", label = TRUE, label.size = 1, 
order=c("Patient3","Other"),
repel = TRUE,pt.size = 0.05) + 
ggtitle("WNN") +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
print(p)
dev.off()


combined$Patient4=rep("Other",length(combined$Label))
combined$Patient4[which(combined$Label=="Patient4")]="Patient4"
Idents(combined)="Patient4"
pdf("Fig1E_Pat4.pdf",width=6,height=5)
p=DimPlot(combined, reduction = "wnn.umap", 
cols=c(alpha(pal_other[1],0.5),alpha(pal_pat[3],0.8)),
group.by = "Patient4", label = TRUE, label.size = 1, 
order=c("Patient4","Other"),
repel = TRUE,pt.size = 0.05) + 
ggtitle("WNN") +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
print(p)
dev.off()


###Normal
combined$Normal=rep("Other",length(combined$Label))
combined$Normal[which(combined$AE=="Normal")]="Normal"
Idents(combined)="Normal"
pdf("Fig1E_Normal_all.pdf",width=6,height=5)
p=DimPlot(combined, reduction = "wnn.umap", 
cols=c(alpha(pal_other[1],0.8),alpha(pal_Nor[1],0.6)),
group.by = "Normal", label = TRUE, label.size = 1, 
order=c("Normal","Other"),
repel = TRUE,pt.size = 0.05) + 
ggtitle("WNN") +
labs(x = "UMAP1", y = "UMAP2") +     
theme(axis.text.y = element_blank(),          
axis.ticks.y = element_blank(),         
axis.text.x = element_blank(),         
axis.ticks.x = element_blank())+    
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))  #加边框
print(p)
dev.off()



6.Fig 1F (cell’s Propotion of each individual_
####Figure 1
library(Seurat)
library(ggplot2)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig1/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#A3D8F4","#8CC0DE","#A1CCD1","#64C9CF","#01A9B4","#176B87")


###Figure 1F,x lab celltype,ylab cell ratio
colnames(combined@meta.data)
cellnumber=apply(table(combined@meta.data[,c(48,60)]),2,sum)
data=table(combined@meta.data[,c(60,48)])/apply(table(combined@meta.data[,c(48,60)]),2,sum)
data1=table(combined@meta.data[,48])/length(combined@meta.data[,48])

###平衡每个样本的细胞数差异性：
data2=c()
for(i in 1:9){
res=data[,i]/data1[i]
data2=cbind(data2,res)
}
colnames(data2)=colnames(data)
rownames(data2)=rownames(data)
data2=as.data.frame(data2)
data2$cell=rownames(data2)
data2=reshape2::melt(data2)
colnames(data2)=c("celltype","Sample","Freq")

data3=table(combined@meta.data[,48])/length(combined@meta.data[,48])
data3=as.data.frame(data3)
data3=cbind(rep("All_Celltypes",length(data3[,1])),data3)
colnames(data3)=c("celltype","Sample","Freq")
data4=rbind(data2,data3)

###把每一个celltype中的比例再缩小到总和为1
cell=unique(data4[,1])
for(n in 1:length(unique(data4[,1]))){
 data4[data4$celltype==cell[n],3]=data4[data4$celltype==cell[n],3]/sum(data4[data4$celltype==cell[n],3])
}

#pal <- paletteer::paletteer_d("palettetown::walrein")[c(9:12,1:3,5:6,8)]
pdf("Fig1F_barplot_cellratio_change.pdf",width=8,height=6)
ggplot(data4,
aes(x=factor(celltype,
             levels=c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
                      "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK",
                       "CM","ncM","cDC","pDC","All_Celltypes")),
    y=as.numeric(Freq),
    fill=factor(Sample,levels=c("Patient1","Patient3","Patient4",
                               "Normal1","Normal2","Normal3","Normal4","Normal5","Normal6"))))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=c(pal_pat[1:3],pal_Nor),name="Samples")+
theme_bw()+
labs(x="Cell Type",y="Proportion")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()