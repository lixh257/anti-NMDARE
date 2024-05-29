
Fig2 Bcell Analysis
###
library(FigR)
library(SimBu)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(plyranges)
library(networkD3)
library(htmltools)
library(ggrepel)
library(Seurat)
library(Signac)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
combined_B=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds") 

1. Fig 2A

pdf("Fig2A_Dimplot_Bmn_celltype.pdf",height=5,width=6)
DimPlot(combined_B, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()


pal <- c("#756AB6","#AC87C5","#E0AED0")
BCell=combined_B
sample_table <- as.data.frame(table(BCell@meta.data$AE,BCell@meta.data$celltype2))
names(sample_table) <- c("Samples","celltype","CellNumber")
sample_table$Samples <- factor(sample_table$Samples, level=c("Normal","Patient"))
sample_table$celltype <- factor(sample_table$celltype, level=c("B naive","B memory","B plasma"))

pdf("Percentage-celltype-in-sample_group.pdf")                     
plot_samplebar<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_samplebar
dev.off()

colnames(BCell@meta.data)
Bp_pro=table(BCell@meta.data[,c(52,60)])[,3]/table(BCell$celltype2)[3]
data=as.data.frame(Bp_pro)
data$Label=c("Normal","Patient")
library(ggplot2)
pdf("Fig2B_Bplasma_propotion.pdf",width=8,height=6)
ggplot(data,
aes(x=factor(Label,levels=c("Normal","Patient")),
y=Bp_pro,
fill=factor(Label,levels=c("Normal","Patient"))))+
geom_bar(stat="identity", width=.8, position = "dodge")+
scale_fill_manual(values=pal[1:2],name="Cell")+
theme_bw()+
labs(x="Label",y="The proportion of each Group")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()


2. Fig 2B (propotion)

colnames(BCell@meta.data)
table(BCell@meta.data[,c(48,60)])
data=table(BCell@meta.data[,c(48,60)])/apply(table(BCell@meta.data[,c(48,60)]),1,sum)
data=as.data.frame(data)
library(ggplot2)
pal=paletteer::paletteer_d("MetBrewer::Benedictus")[c(13,12,11,10,9,8,6,5,4,3,2)]
pdf("Fig2B_Bcell_propotion.pdf",width=8,height=6)
ggplot(data,
aes(x=factor(Label,levels=c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Patient1","Patient3","Patient4")),
y=Freq,
fill=factor(celltype2,levels=c("B naive","B memory","B plasma"))))+
geom_bar(stat="identity", width=.8, position = "dodge")+
scale_fill_manual(values=pal[1:3],name="Cell")+
theme_bw()+
labs(x="Samples",y="The proportion of each Sample")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()


3.Fig 2C
###计算不同celltype-Sample的相关性
library(Signac)
library(Seurat)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds") 
###计算每个细胞之间的相关性
##提取差异较大的2000个基因作为计算相关性的feature
DefaultAssay(BCell)="RNA"
exp=GetAssayData(BCell,slot="data")
features=names(tail(sort(apply(exp, 1, sd)),2000))
exp=exp[which(row.names(exp)%in% features),]
exp=as.data.frame(exp)
cor.exp <- cor(exp, method= "spearman")

#change  cell order
nn1=which(BCell$celltype2=="B naive" & BCell$Label=="Normal1")
nn2=which(BCell$celltype2=="B naive" & BCell$Label=="Normal2")
nn3=which(BCell$celltype2=="B naive" & BCell$Label=="Normal3")
nn4=which(BCell$celltype2=="B naive" & BCell$Label=="Normal4")
nn5=which(BCell$celltype2=="B naive" & BCell$Label=="Normal5")
nn6=which(BCell$celltype2=="B naive" & BCell$Label=="Normal6")
np1=which(BCell$celltype2=="B naive" & BCell$Label=="Patient1")
np2=which(BCell$celltype2=="B naive" & BCell$Label=="Patient3")
np3=which(BCell$celltype2=="B naive" & BCell$Label=="Patient4")

mn1=which(BCell$celltype2=="B memory" & BCell$Label=="Normal1")
mn2=which(BCell$celltype2=="B memory" & BCell$Label=="Normal2")
mn3=which(BCell$celltype2=="B memory" & BCell$Label=="Normal3")
mn4=which(BCell$celltype2=="B memory" & BCell$Label=="Normal4")
mn5=which(BCell$celltype2=="B memory" & BCell$Label=="Normal5")
mn6=which(BCell$celltype2=="B memory" & BCell$Label=="Normal6")
mp1=which(BCell$celltype2=="B memory" & BCell$Label=="Patient1")
mp2=which(BCell$celltype2=="B memory" & BCell$Label=="Patient3")
mp3=which(BCell$celltype2=="B memory" & BCell$Label=="Patient4")

pn1=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal1")
pn2=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal2")
pn3=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal3")
pn4=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal4")
pn5=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal5")
pn6=which(BCell$celltype2=="B plasma" & BCell$Label=="Normal6")
pp1=which(BCell$celltype2=="B plasma" & BCell$Label=="Patient1")
pp2=which(BCell$celltype2=="B plasma" & BCell$Label=="Patient3")
pp3=which(BCell$celltype2=="B plasma" & BCell$Label=="Patient4")

cor.exp1 <- cor.exp[,c(nn1,nn2,nn3,nn4,nn5,nn6,np1,np2,np3,
                       mn1,mn2,mn3,mn4,mn5,mn6,mp1,mp2,mp3,
                       pn1,pn2,pn3,pn4,pn5,pn6,pp1,pp2,pp3)] 
cor.exp1 <- cor.exp1[c(nn1,nn2,nn3,nn4,nn5,nn6,np1,np2,np3,
                       mn1,mn2,mn3,mn4,mn5,mn6,mp1,mp2,mp3,
                       pn1,pn2,pn3,pn4,pn5,pn6,pp1,pp2,pp3),]
cor.exp1=as.data.frame(cor.exp1)
#colume and row annotation #remove normal4 normal5
annotation_col = data.frame(
  Label = c(rep("Normal1",length(nn1)),rep("Normal2",length(nn2)),rep("Normal3",length(nn3)),rep("Normal4",length(nn4)),
            rep("Normal5",length(nn5)),rep("Normal6",length(nn6)),rep("Patient1",length(np1)),rep("Patient3",length(np2)),
            rep("Patient4",length(np3)),
            rep("Normal1",length(mn1)),rep("Normal2",length(mn2)),rep("Normal3",length(mn3)),rep("Normal4",length(mn4)),
            rep("Normal5",length(mn5)),rep("Normal6",length(mn6)),rep("Patient1",length(mp1)),rep("Patient3",length(mp2)),
            rep("Patient4",length(mp3)),
            rep("Normal1",length(pn1)),rep("Normal2",length(pn2)),rep("Normal3",length(pn3)),rep("Normal4",length(pn4)),
            rep("Normal5",length(pn5)),rep("Normal6",length(pn6)),rep("Patient1",length(pp1)),rep("Patient3",length(pp2)),
            rep("Patient4",length(pp3))), 
  celltype = c(rep("B naive",length(which(BCell$celltype2=="B naive"))),
               rep("B memory",length(which(BCell$celltype2=="B memory"))),
               rep("B plasma",,length(which(BCell$celltype2=="B plasma")))) 
  )
row.names(annotation_col) <- colnames(cor.exp1)
annotation_row = annotation_col
row.names(annotation_row) <- rownames(cor.exp1)
# Define custom color palettes for the "Label" and "celltype" annotations in annotation_col
label_colors <- c("B naive" = paletteer::paletteer_d("ggsci::category20_d3")[10],
                  "B memory" = paletteer::paletteer_d("ggsci::category20_d3")[1],
                  "B plasma" = paletteer::paletteer_d("ggsci::category20_d3")[20]
                  )
celltype_colors <- c("Normal1" = "#539092", 
                     "Normal2" = "#8BCFCC", 
                     "Normal3" = "#51ADCF",
                     "Normal4" = "#0278AE",
                     "Normal5" = "#01A9B4",
                     "Normal6" = "#A3D8F4",
                     "Patient1" = "#756AB6",
                     "Patient3" = "#AC87C5",
                     "Patient4" = "#E0AED0"
                    )
# Define custom color palettes for the "Label" and "celltype" annotations in annotation_row
row_label_colors  <- label_colors
row_celltype_colors <- celltype_colors
#change cluster_cols = FALSE, cluster_rows = FALSE to TURE
library(pheatmap)
pdf("Fig2C_Correlation-heatmap-among-celltype-label-cluster-1.pdf",width=8,height=7)
png("Fig2C_Correlation-heatmap-among-celltype-label-cluster-1.png",width=1000,height=1000)
pheatmap(cor.exp1,rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         cluster_cols = FALSE, cluster_rows = FALSE,
         treeheight_row = 5,treeheight_col = 5,
         breaks=(c(seq(0, 0.25, length.out = 7),seq(0.25001, 0.5, length.out = 3))),     
         show_rownames = F,show_colnames = F,       
         border_color = NA,   scale = "none",      
         annotation_colors = list(celltype = label_colors, Label = celltype_colors),
         # Apply custom color palettes for annotation_col
         annotation_colors_row = list(celltype = row_label_colors, Label = row_celltype_colors)
         )
dev.off()



4-5.Fig 2D-E
##################################################################################
#按照celltype添加peak 和chromvar matrix
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds") 
cell=c("B naive","B memory","B plasma")
Idents(BCell)="celltype2"
DefaultAssay(BCell) <- 'ATAC'
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# call peaks for each celltype in myeloid using MACS2
peaks <- CallPeaks(BCell, group.by = "celltype2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks1 <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks1 <- subsetByOverlaps(x = peaks1, ranges = blacklist_hg38_unified, invert = TRUE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(BCell),
  features = peaks1,
  cells = colnames(BCell)
)
fragpath <- Fragments(BCell)
# create a new assay using the MACS2 peak set and add it to the Seurat object
BCell[["Peaks1"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(BCell) <- "Peaks1"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
)
motif.mat <- CreateMotifMatrix(
  features = StringToGRanges(rownames(BCell)),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
# add motif information
BCell <- AddMotifs(
  object = BCell,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
BCell <- RegionStats(
  object = BCell,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)
#RunChromVAR to caculate the motif activities
DefaultAssay(BCell) <- "Peaks1"
# create a peak x motif matrix, where each entry is 1 if the peak contains the motif, otherwise 0
motif.mat <- CreateMotifMatrix(
  features = StringToGRanges(rownames(BCell)),
  pwm = pfm,
  genome = 'hg38'
)
BCell <- RunChromVAR(
  object = BCell,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.mat,
  new.assay.name = "chromvar"
)
saveRDS(BCell,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")


# 找细胞类型之间的差异性
#4.1 RNA
library(SeuratWrappers)
library(SeuratObject)
library(ggplot2)
library(devtools)
library(ggpubr)
library(dplyr)
library(Seurat)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
DefaultAssay(BCell) <- 'RNA'
Idents(BCell) <- "celltype2"
DEG_celltype <- FindAllMarkers(
  object = BCell,
  min.pct = 0.1,
  logfc.threshold=0,  # Default 0.25
  test.use = 'LR',
  only.pos = TRUE
  )
DEG_celltype1 <- DEG_celltype[DEG_celltype$p_val_adj <= 0.01 &DEG_celltype$avg_log2FC >= 0.5,]
write.table(DEG_celltype1,"DEG_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",quote=F)


##用FigR找上面那些DEG附近的高相关的peak
library(FigR)
library(SimBu)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(plyranges)
library(networkD3)
library(htmltools)
library(ggrepel)
library(Seurat)
library(Signac)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
res=BCell
DefaultAssay(res) <- "Peaks1"
counts_1=res@assays$Peaks1@data
res <- RunTFIDF(res,assay="Peaks1")##runGenePeakcorr()里面默认输入的是没有做标准化的ATAC matrix，然后它函数里面用的是centerCounts()做了简单的行归一化
###我们这里用TFIDF做标准化，然后函数里面的参数normalizeATACmat=FALSE，让函数不用再做标准化了
counts=res@assays$Peaks1@data#data就是runTFIDF的结果
#counts = as.matrix(counts)
counts <- counts[Matrix::rowSums(counts)!=0,]
colData <- DataFrame(res@meta.data)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
AllPeak <- ClosestFeature(res, regions = rownames(counts),annotation = annotation,sep = c(':', '-'))
Granges <- data.frame(peak=AllPeak$query_region)
Granges <- str_split_fixed(Granges$peak, ":", 2) 
Granges <- data.frame(chr=Granges[,1],range=Granges[,2])
Granges2 <- str_split_fixed(Granges$range, "-", 2)
Granges <- data.frame(chr=Granges[,1],start=Granges2[,1],end=Granges2[,2])
#head(Granges)
Granges <- cbind(AllPeak,Granges)
Granges$start <- as.numeric(Granges$start)
Granges$end <- as.numeric(Granges$end)
peaks_gr <- as_granges(Granges,seqnames = chr) 
ATAC.SE <- SummarizedExperiment(assays = list(counts=counts), rowRanges = peaks_gr,colData = colData)
##RNA 
DefaultAssay(res) <- 'RNA' 
res <- SCTransform(res, vst.flavor = "v2",method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE) 
DefaultAssay(res) <- 'SCT' 
rnaMat=res@assays$SCT@data###runGenePeakcorr()输入的RNA matrix需要是normalized,这里我们不用RNA@data(normalization)的，用SCT normal 的矩阵
rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]

#Peak-gene association by gene's TSS window and peak accessibility(mean-centered) 
cisCor <- runGenePeakcorr(ATAC.se = ATAC.SE,
                           RNAmat = rnaMat,
                           genome = "hg38",
                           nCores = 16, 
                           normalizeATACmat=FALSE,
                           windowPadSize=200000,
                           p.cut=NULL,keepMultiMappingPeaks=TRUE,
                           n_bg = 100)
# Filter peak-gene correlations by p-value                    
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.01)
write.table(cisCor.filt,"Bcell_all_cisCor_filter_200kb.txt",sep="\t",quote=F)


###intersect the DORC and DEG
DEG=read.table("DEG_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",head=T)
DEG$Gene=DEG$gene
DORC=read.table("Bcell_all_cisCor_filter_200kb.txt",sep="\t",head=T)
DEG_DORC=unique(intersect(DEG$gene,DORC$Gene))
res1=c()
for(i in 1:length(DEG_DORC)){
  print(i)
  res=c()
  fc_res=DEG[which(DEG$gene==DEG_DORC[i]),c(2,5,6,7)]
  cis_res=DORC[which(DORC$Gene==DEG_DORC[i]),2:5]
  for(n in 1:length(fc_res[,1])){
  res=cbind(cis_res,fc_res[n,])
  res1=rbind(res1,res)
}
}
write.table(res1,"DEG_DORC_Cis_peak_res.txt",sep="\t",quote=F)


####重新定义这些交集的DEG为最终的差异基因集合，再进行画图和gene-peak link分析
DEG=read.table("DEG_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",head=T)
DEG_peak=read.table("DEG_DORC_Cis_peak_res.txt",sep="\t",head=T)
#DEP=read.table("DEP_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",head=T)
DEG_peak$peaks=gsub(":","-",DEG_peak$PeakRanges)
gene_cell=unique(DEG_peak[,c(2,7)])
gene=c(gene_cell[gene_cell[,2]=="B naive",1],gene_cell[gene_cell[,2]=="B memory",1],gene_cell[gene_cell[,2]=="B plasma",1])
#length(unique(DEG_peak$Gene))#184
#length(unique(DEG_peak$peaks))#389个
library(SeuratWrappers)
library(SeuratObject)
library(ggplot2)
library(devtools)
library(ggpubr)
library(dplyr)
library(Seurat)
library(Rmagic)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
Idents(BCell) <- "celltype2"
DefaultAssay(BCell) <- 'RNA'
data_rna=GetAssayData(BCell,slot="data")
range1 <- 1:2336 #Naive B
range2 <- 2337:4666 #Memory B
range3 <- 4676:4963 #Plasma B
data_rna1=cbind(data_rna[,which(BCell$celltype2=="B naive")],
               data_rna[,which(BCell$celltype2=="B memory")],
               data_rna[,which(BCell$celltype2=="B plasma")])
# Set the number of columns to select from each range
num_select <- 100
# Randomly select columns from each range
selected_cols_range1 <- sample(range1, num_select)
selected_cols_range2 <- sample(range2, num_select)
selected_cols_range3 <- sample(range3, num_select)
data_rna2=data_rna1[,c(selected_cols_range1,selected_cols_range2,selected_cols_range3)] 
gene1=intersect(gene,rownames(data_rna2)) 
data_rna3=data_rna2[gene1,]           
library(pheatmap)
library(paletteer)
showname=c("CD27","CD38","MZB1","BCL2L11","JCHAIN","IFNAR2","IL6ST","CDK6","ELL2","ZBTB20","HSPA5",
           "FNDC3B","IRF4","CDK17","PRDM1","ERN1","EIF3A","TCF3","IRF1","IRF9","XBP1")#
lab_row=ifelse(rownames(data_rna3)%in%showname,rownames(data_rna3),"")
pheatmap(data_rna3,labels_row=lab_row,
         color = colorRampPalette(paletteer::paletteer_d("khroma::BuRd"))(100),
         breaks=(c(seq(-2, 0, length.out = 50),seq(0.001, 2, length.out = 50))),
         cluster_rows = FALSE,   cluster_cols = FALSE,
         show_rownames = TRUE,show_colnames = FALSE, #annotation_row = ann_row,
         #annotation_col = column_name,  annotation_colors = annoColor,
         border_color = NA,fontsize = 4,scale = "row",   
         file="pheatmap_DEG_celltype2_2.pdf",###Fig2D
         width=4,height=6)

peak_res=c()
for(g in 1:length(gene)){
  peak_res=rbind(peak_res,DEG_peak[DEG_peak$Gene==gene[g],8:9])
}
peak=peak_res[,2]

#4.2 peaks
library(Rmagic)
Idents(BCell) <- "celltype2"
DefaultAssay(BCell) <- 'Peaks1'
data_peak=GetAssayData(BCell,slot="data")
data_peak_magic=magic(data_peak)
data_peak_magic=data_peak_magic$result
data_peak_magic1=cbind(data_peak_magic[,which(BCell$celltype2=="B naive")],
                 data_peak_magic[,which(BCell$celltype2=="B memory")],
                 data_peak_magic[,which(BCell$celltype2=="B plasma")])
peak1=intersect(peak,rownames(data_peak)) 
data_peak_magic2=data_peak_magic1[peak1,] 
range1 <- 1:2336 #Naive B
range2 <- 2337:4666 #Memory B
range3 <- 4676:4963 #Plasma B
# Set the number of columns to select from each range
num_select <- 100
# Randomly select columns from each range
selected_cols_range1 <- sample(range1, num_select)
selected_cols_range2 <- sample(range2, num_select)
selected_cols_range3 <- sample(range3, num_select)
data_peak2=data_peak_magic2[,c(selected_cols_range1,selected_cols_range2,selected_cols_range3)]           
library(pheatmap)
library(paletteer)
pheatmap(data_peak2,
         color = colorRampPalette(c("#4393C3FF","#92C5DEFF","#D1E5F0FF","white","#F4DE3AFF","#FCB11CFF","#E6B91EFF"))(100),
         #color = colorRampPalette(paletteer::paletteer_d("khroma::BuRd"))(100),
         breaks=(c(seq(-2, 0, length.out = 50),seq(0.001, 2, length.out = 50))),
         cluster_rows = FALSE,   cluster_cols = FALSE,
         show_rownames = FALSE,show_colnames = FALSE, #annotation_row = ann_row,
         #annotation_col = column_name,  annotation_colors = annoColor,
         border_color = NA,fontsize = 4,scale = "row",   
         file="pheatmap_DEP_celltype2_1.pdf",###Fig2E
         width=4,height=6)
#write.table(peak_res,"DEG_cis_peak_res.txt",sep="\t",quote=F)



6.Fig 2F
###Function enrichment
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
##184 gene -metascape-整理过高富集的那些，并且删掉了一些重复的（在不同数据库中相似的term只保留一个）
res=read.table("Bp_184gene_Term_result.txt",sep="\t",head=T)
res$Pvalue= -res$LogP
res$Ratio=res$InTerm_count/res$Term_count
res$Description=factor(res$Description,levels=res$Description)

library(ggplot2)
pdf("Fig2F_Bp_184gene_Founction_enrichment_res.pdf",width=9,height=4)
ggplot(res,aes(x=Pvalue,y=Description,color=Ratio,size=InTerm_count))+geom_point()+
scale_colour_gradientn(colours =paletteer::paletteer_d("RColorBrewer::PuOr")[7:11] )+
theme_bw()
dev.off() 


7-8. Fig2G-H Motif Analysis
###Rank Motif
####Bplasma feature plot
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
DEG_peak=read.table("DEG_DORC_Cis_peak_res.txt",sep="\t",head=T)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(AnnotationDbi)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
## Plasma B markerPeak
PlasmaB_markerPeak <- DEG_peak[DEG_peak$cluster =="B plasma",]
PlasmaB_markerPeak$peaks=gsub(":","-",PlasmaB_markerPeak$PeakRanges)
## make GRanger Format of MarkerPeak list
PlasmaB_Granges <- str_split_fixed(PlasmaB_markerPeak$peaks, "-", 2) 
PlasmaB_Granges <- data.frame(chr=PlasmaB_Granges[,1],range=PlasmaB_Granges[,2])
PlasmaB_Granges2 <- str_split_fixed(PlasmaB_Granges$range, "-", 2)
PlasmaB_Granges <- data.frame(chr=PlasmaB_Granges[,1],start=PlasmaB_Granges2[,1],end=PlasmaB_Granges2[,2])
PlasmaB_markerPeak_site <- cbind(PlasmaB_markerPeak,PlasmaB_Granges)
PlasmaB_markerPeak_site$start <- as.numeric(PlasmaB_markerPeak_site$start)
PlasmaB_markerPeak_site$end <- as.numeric(PlasmaB_markerPeak_site$end)
PlasmaB_markerPeak_site <- makeGRangesFromDataFrame(PlasmaB_markerPeak_site)
markerpeak_gene <- annotatePeak(PlasmaB_markerPeak_site, TxDb=txdb,tssRegion=c(-3000, 3000), 
                                addFlankGeneInfo=TRUE, flankDistance=200000, annoDb="org.Hs.eg.db")
#markerpeak_gene
write.csv(data.frame(markerpeak_gene),"DEG_PlasmaB-up-cisPeak.csv",quote=F)

pdf("PlasmaB-up-markerPeak-distribution.pdf")
plotAnnoBar(markerpeak_gene)
plotDistToTSS(markerpeak_gene,title="Distribution of transcription MarkerPeak loci \n relative to TSS")
dev.off()
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(dplyr)
###subset Bp cell and findmotif
Bp=subset(BCell,celltype2=="B plasma")
DefaultAssay(Bp)="Peaks1"
Bp <- RegionStats(Bp, genome = BSgenome.Hsapiens.UCSC.hg38)
PlasmaB_markerPeak_motif <- FindMotifs(Bp, features = PlasmaB_markerPeak$peaks)
# Rank motifs by p-value and fold enrichment
PlasmaB_markerPeak_motif_rank <- PlasmaB_markerPeak_motif[order(PlasmaB_markerPeak_motif$pvalue, PlasmaB_markerPeak_motif$fold.enrichment), ]
write.table(PlasmaB_markerPeak_motif_rank,"1-PlasmaB_markerPeak_motif_rank.txt",sep="\t",quote=F)
#extract motif name
Motif_TF_name <- PlasmaB_markerPeak_motif_rank[,c(1,8)]
write.table(Motif_TF_name,"Motif_TF_name.txt",sep="\t",quote=F)
##### Specific motif by differential chromvar score #####
DefaultAssay(BCell) <- 'chromvar'
Idents(BCell) <- "celltype2"
Bplasma_chromvarMotif_vsOtherB <- FindMarkers(
  object = BCell,
  ident.1="B plasma",
  ident.2=c("B memory","B naive"),
  min.pct = 0.1,
  logfc.threshold=0,  # Default 0.25
  #only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
Bplasma_chromvarMotif_vsOtherB$motif=rownames(Bplasma_chromvarMotif_vsOtherB)
Motif_res1=merge(PlasmaB_markerPeak_motif_rank,Bplasma_chromvarMotif_vsOtherB,by = "motif")
Motif_res1=Motif_res1[order(Motif_res1$pvalue),]
Motif_res1=Motif_res1[,c(8,1,4,6,9,14,11)]
colnames(Motif_res1)=c("motif.name","motif","enrich_percent.observed",
                       "enrich_fold.enrichment","enrich_padj","chromvar_padj","chromvar_avg_diff")

#add expression of the motif
DefaultAssay(BCell) <- "RNA"
Idents(BCell) <- "celltype2"
aa=GetAssayData(BCell)
ave.exp=AverageExpression(BCell)
exp_mean=ave.exp$RNA
###motif.name
Motif=read.table("Motif_TF_name.txt",sep="\t",head=T)###手动修改这个cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5_Bcell/ /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2
#setdiff(motif$motif.name.1,rownames(exp_mean))
Motif_res1$RNA_log2FC=rep(1,length(Motif_res1[,1]))
Motif_res1$Bp_mean_RNA=rep(0,length(Motif_res1[,1]))
Motif_res1$Bp_RNA_ratio=rep(0,length(Motif_res1[,1]))
for(m in 1:length(Motif_res1[,1])){
  tf_gene=Motif[which(Motif[,2]==Motif_res1[m,2]),4]
  RNA_FC=(exp_mean[tf_gene,3]+0.01)/((exp_mean[tf_gene,1]+exp_mean[tf_gene,2])/2+0.01)
  RNA_FC_m=sum(RNA_FC)/length(RNA_FC)
  Motif_res1$RNA_log2FC[m]=log2(RNA_FC_m)
  RNA_mean=exp_mean[tf_gene,3]
  RNA_mean_m=sum(RNA_mean)/length(RNA_mean)
  Motif_res1$Bp_mean_RNA[m]=RNA_mean_m
  RNA_ratio_m=length(which(aa[tf_gene,which(BCell$celltype2=="B plasma")]>0))/length(which(aa[tf_gene,which(BCell$celltype2=="B plasma")]>=0))
  Motif_res1$Bp_RNA_ratio[m]=RNA_ratio_m
}
rownames(Motif_res1) <- Motif_res1$motif
write.csv(Motif_res1,"DEG_PlasmaB-up-cisPeak_Motif_res.csv",quote=F)

##Ranking Motifs by enrich pvalue
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
res1=read.csv("DEG_PlasmaB-up-cisPeak_Motif_res.csv")
library(ggplot2)
library(ggrepel)
#res=Bplasma_chromvarMotif_vsOtherB2
res1$Rank=1:length(res1[,1])
#res1=res1[1:length(which(res1$enrich_padj<0.01)),]
##mark gene
####
pdf("Fig2G_Motif-ranking_Bplasma.pdf",height=4,width=10)
p1 <- ggplot(res1,
  aes(x=Rank, y=-log10(enrich_padj),color=enrich_fold.enrichment,
      size=enrich_percent.observed)) +
	geom_point(alpha=0.8) + #, size=2
	labs(x= "Motif Rank") +ylab("-log10(pval)") +
	scale_color_gradient(low = "yellow", high = "darkblue")+
	geom_text_repel(data = res1[which(res1$enrich_padj<=0.00001),], #
  aes(label = motif.name),max.overlaps=Inf, vjust = -0.5,color="black",size=2,
  box.padding=unit(0.5,'lines'), point.padding=unit(0.1, 'lines'), 
  segment.color='black', show.legend=FALSE
  )+theme_bw()+
  scale_size(breaks = c( 10,20, 30,40),range = c(0,6))
print(p1)
dev.off()


pdf("Fig2G_Motif-ranking_Bplasma_new.pdf",height=4,width=8)
p1 <- ggplot(res1,
  aes(x=enrich_fold.enrichment, y=-log10(enrich_padj),color=enrich_fold.enrichment,
      size=enrich_percent.observed)) +
	geom_point(alpha=0.8) + #, size=2
	labs(x= "enrichment score") +ylab("-log10(padj)") +
	scale_color_gradient(low = "yellow", high = "darkblue")+
	geom_text_repel(data = res1[which(res1$enrich_padj<=0.05),], #
  aes(label = motif.name),max.overlaps=Inf, vjust = -0.5,color="black",size=2,
  box.padding=unit(0.5,'lines'), point.padding=unit(0.1, 'lines'), 
  segment.color='black', show.legend=FALSE
  )+
  theme_bw()+
  scale_size(breaks = c( 10,20, 30,40),range = c(0,6))
print(p1)
dev.off()



res_mark1=res1[which(res1$enrich_padj<=0.05 &res1$RNA_log2FC>0),]
pdf("Fig2H_Motif-ranking_Bp_FC_new.pdf",height=5,width=9)
p1 <- ggplot(res_mark1,
  aes(x=chromvar_avg_diff, y=Bp_RNA_ratio,
  color=RNA_log2FC,size=-log10(chromvar_padj+ 1e-310))) +
	geom_point(alpha=1) + #, size=2
	labs(x= "TF activity(log2FC)") +ylab("TF expression(Ratio)") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
	geom_text_repel(data = res_mark1, 
  aes(label = motif.name), max.overlaps=Inf,vjust = -0.5,color="black",size=3.5,
  box.padding=unit(0.5,'lines'), point.padding=unit(0.1, 'lines'), 
  segment.color='black', show.legend=FALSE
  )+#xlim(c(-1,1))+ylim(c(-0.5,1))+
  theme_bw()+
  scale_size(breaks = c(10,20,30,50,100,300),range = c(1,10))
print(p1)
dev.off()



9.Fig 2I GRN
###整理H图中的TF的调控基因list
###查看Bindingsite
#用MODDS()计算Binding Score 和 P值以及Binding sites
##把pwm矩阵输出成单独的pfm文件
pwm_all=Bp@assays$Peaks1@motifs@pwm
MotifID=read.table("Motif_TF_name.txt",sep="\t",head=T)
for(t in 1:length(MotifID[,1])){
 pwm1 = pwm_all[t][[1]]
 write.table(pwm1,paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/633_Motifs_pfm/",
                         MotifID[t,1],"-",MotifID[t,2],"_pwm.pfm"),
                         row.names=F,col.names=F,quote=F)
}
##Peak list 输出成fa文件
DEG_peak=read.table("DEG_DORC_Cis_peak_res.txt",sep="\t",head=T)
DEG_peak$peaks=gsub(":","-",DEG_peak$PeakRanges)
## Plasma B markerPeak
PlasmaB_markerPeak <- DEG_peak[DEG_peak$cluster =="B plasma",]## 365peaks
PlasmaB_markerPeak_bed=gsub("-","\t",PlasmaB_markerPeak$peaks,perl=T)
write.table(PlasmaB_markerPeak_bed,"PlasmaB_DORC_link_peaks_res.bed",sep="\t",row.names=F,col.names=F,quote=F)

#提取peak的DNA序列出来
#linux
source activate Pyth3_R4
bedtools getfasta -fi /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-bed /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/PlasmaB_DORC_link_peaks_res.bed \
-fo /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/PlasmaB_DORC_link_peaks_res.fa

###计算Binding
#linux
cd /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/633_Motifs_pfm/
mkdir /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_Bplasma/
for i in *_pwm.pfm
do
echo "$i"
moods-dna.py --m "$i" \
-s /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/PlasmaB_DORC_link_peaks_res.fa \
--p-value 0.00001  \
--outfile /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_Bplasma/${i%.*}_motif_bind_0.00001.txt
done

######
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
DORC_res=read.table("DEG_DORC_Cis_peak_res.txt",sep="\t",head=T)
DORC_res_Bp=DORC_res[DORC_res$cluster=="B plasma",]
motif_res=read.csv("DEG_PlasmaB-up-cisPeak_Motif_res.csv")
###挑出H图中的TF的结果：
TF=c("IRF9","IRF4","IRF1","IRF7","IRF3","IRF2","STAT1::STAT2",
     "TCF12(var.2)","TWIST1","HSF4")
motif_res_significant=motif_res[motif_res$motif.name%in%TF,]
rownames(motif_res_significant)=motif_res_significant$motif.name

Bind_res4=c()
for(t in 1:10){
  TF_ID1=motif_res_significant[TF[t],1]
  TF_name=motif_res_significant[TF[t],2]
  Bind_res=read.table(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_Bplasma/",
                              TF_ID1,"-",TF_name,"_pwm_motif_bind_0.0001.txt"),sep=",")
  Bind_res1=Bind_res[,c(1,3,4,5,6)]
  colnames(Bind_res1)=c("PeakRanges","width","strand","BindingScore","BindingSeq")
   if(max(Bind_res1$BindingScore)>7){
  Bind_res1=Bind_res1[Bind_res1$BindingScore>=7,]
  peak_gene_1=c()
  for(p in 1:length(Bind_res1[,1])){
    peak_gene=DORC_res_Bp[which(DORC_res_Bp$PeakRanges==Bind_res1[p,1]),1:7]
    peak_gene_1=rbind(peak_gene_1,peak_gene)
  }
  Bind_res2=merge(Bind_res1,peak_gene_1,by="PeakRanges")
  Bind_res3=cbind(Bind_res2,motif_res_significant[TF[t],])
  Bind_res4=rbind(Bind_res4,Bind_res3)
}
}
Bind_res4=Bind_res4[,-12]
colnames(Bind_res4)=c("PeakRanges","Bindingwidth","Bindingstrand","BindingScore","BindingSeq","Gene",
                      "Peak-Gene-rObs","Peak-Gene-pvalZ","Gene-avg_log2FC","Gene-p_val_adj","cluster",
                      "motif.name","motifID","enrich_percent.observed","enrich_fold.enrichment","enrich_padj",
                      "chromvar_padj","chromvar_avg_diff","TF_RNA_log2FC")
write.table(Bind_res4,"Fig2I_DEG_PlasmaB-up-cisPeak_Motif_res_Binding_p_0.0001_result_last.txt",
            sep="\t",quote=F,row.names=F)
###Fig2I network of gene-peak-tf 导入上述表格到CytoScape中
##然后在AI里面调整细节



10.Fig 2J
#######上游的触发通路是什么呢？
####上述IRFs多来源于IFN pathway下游激活，且已知其他自身免疫性疾病和IFN相关
##验证IFN pathway中基因在B细胞中的表达
##IFN gene in Bplasma cell
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(AnnotationDbi)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
###先检验IFN pathway
IFNR=c("IFNAR1","IFNAR2","TYK2","JAK1","JAK2","IL10RB","IFNLR1","STAT1","STAT2","IRF9","IFNGR1","IFNGR2")
DefaultAssay(BCell) <- 'RNA' 
RNA_data=GetAssayData(BCell,slot="data")
RNA_data1=RNA_data[IFNR,]
RES=c()
for(i in 1:length(IFNR)){
  Bp_exp=RNA_data1[IFNR[i],which(BCell$celltype2=="B plasma")]
  Bnm_exp=RNA_data1[IFNR[i],which(BCell$celltype2!="B plasma")]
  Test=wilcox.test(Bp_exp,Bnm_exp)
  Test_p=Test$p.value
  FC=mean(Bp_exp)/mean(Bnm_exp)
  res=cbind(IFNR[i],Test_p,FC)
  RES=rbind(RES,res)
}
##pointplot
pdf("Fig2J_IFNR_gene_Bplasma_exp_res.pdf",width=5,height=5)
DotPlot(BCell, features = IFNR,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,10, 20, 30, 50,75),range = c(-1,7))+
labs(y="Cell Type",x="IFNR Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma")
)
dev.off()


11.Fig 2K 对于IFN pathway，画模式图展示Bplasma中的调控关系





























##vlnplot
DefaultAssay(BCell)="RNA"
Idents(BCell)="celltype2"
levels(BCell) <- c("B naive","B memory","B plasma")
pdf(paste0("FigS_Bcell_IFN_Pathway_gene_RNA.pdf"),width=10,height=5)
for(i in 1:length(IFNR)){
p=VlnPlot(BCell, features =IFNR[i], 
         pt.size = 0.1, 
        combine = T, #split.plot = T,
        cols =paletteer::paletteer_d("Polychrome::palette36") 
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("CellType") + ylab("Gene Expression") + 
       #stat_compare_means( method = "wilcox.test") +
       #stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()

Bp=subset(BCell,celltype2=="B plasma")
Idents(Bp)="Label"
levels(Bp)=c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Patient1","Patient3","Patient4")

pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#539092","#8BCFCC","#51ADCF","#0278AE","#01A9B4","#A3D8F4")

pdf(paste0("FigS_Bcell_IFN_Pathway_gene_RNA_PvsN.pdf"),width=6,height=3)
for(i in 1:length(IFNR)){
p=VlnPlot(Bp, features = IFNR[i], 
         pt.size = 0.1, 
        combine = T
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("CellType") + ylab("Gene Expression") + 
       scale_fill_manual(values=c(pal_Nor,pal_pat),name="Samples")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()
























setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
DefaultAssay(BCell)="RNA"
aa=GetAssayData(BCell)
allgene=rownames(aa)

#IFNI pathway:
IFN1=c("IFNA1","IFNB1","TYK2","JAK1","JAK2")
IFN3=c("IL10RB","IFNLR1")
TNF=c("TNFRSF9","TNFSF9","TNFSF13","TNFSF13B","TNFRSF13B","TNFRSF13C","TNFRSF17",
      "CD70","TNFRSF7","TNFRSF8","TNFSF8","CD40","CD40LG",)
TNF[which(TNF%in%allgene)]

DEG=read.table("DEG_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",head=T)
DEG[which(DEG$gene%in%TNF),]
###上面这些通路相关蛋白并没有显著表达变化

setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(AnnotationDbi)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
###先检验IFN pathway
gene=c(
  #IFN
  "IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA12P","IFNA13",#not found:"IFNA11P"
  "IFNA14","IFNA16","IFNA17","IFNA20","IFNA21","IFNA22P",
  "IFNB1","IFNG","IFNL1","IFNL2","IFNL3","IFNNP1","IFNWP2","IFNWP4","IFNWP5","IFNWP9","IFNWP15",
  "IFNWP18","IFNWP19","IFNW1","IL6",
  "TYK2","JAK1","JAK2","IL10RB","IFNLR1","STAT1","STAT2","STAT3","STAT4","STAT5",
  #IFN receptor
  "IFNAR1","IFNAR2",
       #other:
       "CXCL13",
       "TNFRSF9","TNFSF9","TNFSF13","TNFSF13B","TNFRSF13B","TNFRSF13C","TNFRSF17",
      "CD70","TNFRSF8","TNFSF8","CD40","CD40LG",#"TNFRSF7",
      "CXCL9","CXCL10","CXCL11","IRF9","MYD88","IRF4","IRF1","IRF2","TCF12")
DefaultAssay(BCell) <- 'RNA' 
SCT=GetAssayData(BCell,slot="data")
Gene=gene[which(gene%in%rownames(SCT))]
gene=Gene

DefaultAssay(BCell)="Peaks1"
BCell <- RegionStats(BCell, genome = BSgenome.Hsapiens.UCSC.hg38)
BCell <- LinkPeaks(
  object = BCell,
  peak.assay = "Peaks1",
  expression.assay = "RNA",
  genes.use = gene,
  pvalue_cutoff = 0.01,
  score_cutoff = 0
  )
Idents(BCell) <- "celltype2"
levels(BCell) <- c("B naive","B memory","B plasma")
DefaultAssay(BCell)="Peaks1"
#Idents(BCell)="cellAE"
pdf("Fig2L_Pathway_gene_track_200kb_RNA.pdf",height=6,width=10)
for(i in 1:length(gene)){
p=CoveragePlot(
   object = BCell,
   region = gene[i], 
   #ranges = Peaks1,
   features = gene[i],
   expression.assay = "RNA",
   extend.upstream = 200000,
   extend.downstream = 200000,
   #ymax = 20,##change normalized signal range
   links=T,
   window = 200
 )
 print(p)
}
 dev.off()







###################################################################################################################
####整理这些通路中基因的表达
#FigA:



       

####IFN gene
IFN=c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA12P","IFNA13",#not found:"IFNA11P"
  "IFNA14","IFNA16","IFNA17","IFNA20","IFNA21","IFNA22P",
  "IFNB1","IFNG","IFNL1","IFNL2","IFNL3","IFNNP1","IFNWP2","IFNWP4","IFNWP5","IFNWP9","IFNWP15",
  "IFNWP18","IFNWP19","IFNW1")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
DefaultAssay(combined) <- 'RNA' 
RNA_data=GetAssayData(combined,slot="data")
IFN=IFN[IFN%in%rownames(RNA_data)]
RNA_data1=RNA_data[IFN,]

cell=unique(combined$celltype2)
RES=c()
for(i in 1:length(IFN)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[IFN[i],which(combined$celltype2==cell[n]&combined$AE=="Patient")]
  exp_Nor=RNA_data1[IFN[i],which(combined$celltype2==cell[n]&combined$AE=="Normal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  res=cbind(IFN[i],cell[n],Test_p,FC)
  RES=rbind(RES,res)
}
}

##pointplot
Idents(combined)="celltype2"
pdf("FigSBcell_B_IFN_gene_cell_exp_res_1.pdf",width=8,height=10)
DotPlot(combined, features = IFN,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,10, 20, 30, 50,75),range = c(-1,7))+
labs(y="Cell Type",x="IFNR Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits =  c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)
dev.off()




#####TLR gene 
####上游通路是什么？
####上游的通路应该是TLR，那么TLR在哪些细胞中表达较高呢？
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")






combined_M=subset(combined,celltype2%in%c("CM","ncM","cDC"))
combined_M$cellAE=paste0(combined_M$celltype2,"_",combined_M$AE)
Idents(combined_M)="cellAE"
DefaultAssay(combined_M)="RNA"
RNA_data=GetAssayData(combined_M,slot="data")
RNA_data1=RNA_data[gene,]
cell=c("CM","ncM","cDC")
combined$cellAE=paste0(combined$AE,"_",combined$celltype2)
Idents(combined)="cellAE"
pdf(paste0("Fig2L_Pathway_gene_RNA_Pat_vs_Normal.pdf"),width=10,height=5)
for(i in 1:length(gene)){
p=VlnPlot(combined, features = gene[i], 
         pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols =paletteer::paletteer_d("Polychrome::palette36") 
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("CellType") + ylab("Gene Expression") + 
       #stat_compare_means( method = "wilcox.test") +
       #stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()

gene1=c("TLR2","TLR4","TLR1","TLR6","TIRAP","MYD88","IRAK4","IRAK1","IRAK2","TRAF6","TAB2","TAB3","MAP3K7","MAP2K3","MAP2K6",
        "MAP2K4","MAP2K7","MAPK14","MAPK8","NFKB1","NFKB2","REL","RELA","CREB1","JUN",
        "CHUK")#"MKK3-MAP2K3","MKK6",JNK1-MAPK8,"IKBKB"

pdf(paste0("FigSBcell_C_TLR_gene_RNA_Pat_vs_Normal_Myeloid.pdf"),width=6,height=3)
for(i in 1:length(gene1)){
p=VlnPlot(combined_M, features = gene1[i], 
         pt.size = 0.1, 
        combine = T, #split.plot = T,
        cols =paletteer::paletteer_d("Polychrome::palette36") 
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("CellType") + ylab("Gene Expression") + 
       #stat_compare_means( method = "wilcox.test") +
       #stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()


###加Test
gene1=c("TLR2","TLR4","TLR1","TLR6","TIRAP","MYD88","IRAK4","IRAK1","IRAK2","TRAF6","TAB2","TAB3","MAP3K7","MAP2K3","MAP2K6",
        "MAP2K4","MAP2K7","MAPK14","MAPK8","NFKB1","NFKB2","REL","RELA","CREB1","JUN",
        "CHUK")
DefaultAssay(combined_M)="RNA"
RNA_data=GetAssayData(combined_M,slot="data")
RNA_data1=RNA_data[gene1,]
cell=c("CM","ncM","cDC")
RES=c()
for(i in 1:length(gene1)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[gene1[i],which(combined_M$celltype2==cell[n]&combined_M$AE=="Patient")]
  exp_Nor=RNA_data1[gene1[i],which(combined_M$celltype2==cell[n]&combined_M$AE=="Normal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  res=cbind(gene1[i],cell[n],Test_p,FC)
  RES=rbind(RES,res)
}
}
RES1=RES[which(as.numeric(RES[,3])<0.05),]
gene_TLR=unique(RES[,1])

DefaultAssay(combined_M)="RNA"
RNA_data=GetAssayData(combined_M,slot="data")
Idents(combined_M)="cellAE"
levels(combined_M)=c("CM_Normal","ncM_Normal","cDC_Normal","CM_Patient","ncM_Patient","cDC_Patient")
pdf("FigSBcell_C_TLR_gene_RNA_Pat_vs_Normal_Myeloid_1.pdf",width=6,height=8)
DotPlot(combined_M, features = gene_TLR,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,10, 20, 30, 50,75),range = c(-1,7))+
labs(y="Cell Type",x="IFNR Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits =  c("CM_Normal","ncM_Normal","cDC_Normal","CM_Patient","ncM_Patient","cDC_Patient")
)
dev.off()

