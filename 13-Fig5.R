Fig5

1.Fig 5A

library(Seurat)
library(ggplot2)
library(dittoSeq)
library(MySeuratWrappers)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=table(combined_M$celltype2)
pdf("Fig5A_Dimplot_M_celltype.pdf",width=8,height=5)
DimPlot(combined_M, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black")
dev.off()

###右上角
##barplot of every individual cell's propotion
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5/")
colnames(combined_M@meta.data)
sample_table <- as.data.frame(table(combined_M@meta.data$AE,combined_M@meta.data$celltype2))
names(sample_table) <- c("Samples","celltype","CellNumber")
sample_table$Samples <- factor(sample_table$Samples, level=c("Normal","Patient"))
sample_table$celltype <- factor(sample_table$celltype, level=c("CM","ncM","cDC","pDC"))

pal=c("#2CA02CFF" ,"#98DF8AFF","#DBDB8DFF","#E377C2FF","#AEC7E8FF","#F7B6D2FF","#9467BDFF","#7F7F7FFF","#176B87")
pdf("Fig5A_Percentage-celltype-in-sample_group.pdf")                     
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


2.Fig 5B DEG expression pheatmap
##DEG gene
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/")
##1 Find DEG (Patient vs Normal)     
library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(Seurat) 
library(tidyverse)
library(magrittr)
library(scater) 
library(cowplot) 
library(Matrix.utils) 
library(edgeR) 
library(ggrepel)
library(dplyr) 
library(magrittr) 
library(Matrix) 
library(purrr) 
library(reshape2) 
library(S4Vectors) 
library(tibble) 
library(SingleCellExperiment) 
library(pheatmap) 
library(apeglm) 
library(png) 
library(DESeq2) 
library(RColorBrewer)
library(harmony)
library(scuttle)
library(scran)
library(glmGamPoi)
#findmarker( test.use  =  DESeq2  )#这个默认的estimateSizeFactors方法用的是默认的type=ratio,不符合低离散度的数据（0值很多，在某个group里面全是0会报错）
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
DefaultAssay(combined) <- 'RNA' 
Idents(combined) <- 'AE' 
cell=unique(combined$celltype2)
#cell=c("cM","C5_M","ncM","cDC","pDC")
for (i in 1:length(cell)){
res=subset(combined,celltype2==cell[i])
res <- SCTransform(res, vst.flavor = "v2",method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE) 
DefaultAssay(res) <- 'SCT' 
counts=GetAssayData(object = res, slot = "count")
counts=counts[which(rowSums(counts>0)>length(counts[1,])/10),]

metadata <- res@meta.data # Set up metadata as desired for aggregation and DE analysis 
metadata$groups_id <- factor(res$AE) 
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata) 
#DESeqDataSetFromMatrix
all(rownames(metadata) == colnames(counts))  
#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(counts,
colData = metadata,                               
design = ~groups_id)
dds <- estimateSizeFactors(dds, type = 'poscounts',locfunc = stats::median)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
#dds <- nbinomWaldTest(dds)
dds <- nbinomLRT(dds,reduced= ~1)
contrast <- c("groups_id", "Patient", "Normal")  
# resultsNames(dds) 
DESeq_res <- results(dds,contrast=contrast,name="groups_id_Patient_vs_Normal",alpha = 0.99, pAdjustMethod="BH")##??????ֵ????ȫ?????
DESeq_res = na.omit(DESeq_res)###?Ϊcountֵ?0?????NA?????
gene=rownames(sce)
res_tbl <- DESeq_res %>%         
data.frame() %>%         
rownames_to_column(var="gene") %>%         
as_tibble()  
res_table_thres=res_tbl
res_table_thres[,8]="Stable"
res_table_thres[which((res_table_thres[,3] > 0.5)&(res_table_thres[,7] < 0.01)),8]="UP"
res_table_thres[which((res_table_thres[,3] < -0.5)&(res_table_thres[,7] < 0.01)),8]="Down" 
colnames(res_table_thres)[8]="threshold"
write.csv(res_table_thres, paste0("DEGs_",cell[i],"_res.csv"), sep="\t",quote =FALSE,row.names = FALSE)


###############################1.MAplot
logcounts1=log2(counts+1)
mean_count=apply(logcounts1,1,mean)
MA_res=matrix(0,ncol=5,nrow=length(rownames(DESeq_res)))
MA_res[,1]=rownames(DESeq_res)
MA_res[,2]=mean_count[rownames(DESeq_res)]
MA_res[,3]=DESeq_res[,2]
MA_res[,4]=DESeq_res[,6]
MA_res[,5]="Stable"
MA_res[which(DESeq_res[,6]<0.01&abs(DESeq_res[,2])>0.25),5]="DEGs"

colnames(MA_res)=c("gene","Mean_counts","log2FC","Padj","State")
MA_res=as.data.frame(MA_res)
MA_res[,2]=as.numeric(MA_res[,2])
MA_res[,3]=as.numeric(MA_res[,3])
MA_res[,4]=as.numeric(MA_res[,4])
p=ggplot(MA_res,
aes(x = as.numeric(MA_res$Mean_counts), y = as.numeric(MA_res$log2FC), colour=MA_res$State)) +     
geom_point(alpha=0.4, size=2.0) +
scale_color_manual(values=c("red","grey"))+
geom_hline(yintercept = c(-0.005,0.005),lty=4,col="black",lwd=0.8) +     
ggtitle("MAplot of Deseq") +  
xlab("mean(log2(raw_count+1))") +      
ylab("log2FC") +  
theme_bw()
MA_res$label=ifelse(as.numeric(MA_res$Mean_counts) > 0.5 & MA_res$State=="DEGs" ,MA_res$gene,"")
p+geom_text_repel(data = MA_res, aes(x = as.numeric(MA_res$Mean_counts), y = as.numeric(MA_res$log2FC), label=MA_res$label),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)            
ggsave(paste0(cell[i],"_DEGs_Deseq_MAplot.pdf"))

############################################## 2.Volcano plot 
p=ggplot(res_table_thres,
aes(x = res_table_thres$log2FoldChange, y = -log10(res_table_thres$padj), colour = threshold,label=res_table_thres$gene)) +     
geom_point(alpha=0.4, size=2.0) +
scale_color_manual(values=c("blue", "grey","red"))+
#xlim(c(-4, 4)) +
geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +     
ggtitle("Volcano plot of AE vs Normal") +  
xlab("log2 fold change") +      
ylab("-log10 padj-value") +  
theme_bw()+   
#scale_y_continuous(limits = c(0,60)) +     
theme(legend.position = "right",           
plot.title = element_text(size = rel(1.5), hjust = 0.5),           
axis.title = element_text(size = rel(1.25))) 

res_table_thres$label=ifelse(res_table_thres$padj < 1e-2 & abs(res_table_thres$log2FoldChange) >= 0.25,res_table_thres$gene,"")

p+geom_text_repel(data = res_table_thres, aes(x = log2FoldChange, 
                                   y = -log10(padj), 
                                   label = label),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)                                                 
               
ggsave(paste0(cell[i],"_DEGs_Volcanoplot_log2FC_0.25.pdf"))
print(cell[i])
}


##2 subset cell's rds and call peak by AE
#####Myeloid
library(FigR)
library(SimBu)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(plyranges)
library(networkD3)
library(htmltools)
library(ggrepel)
library(chromVARmotifs)
library(dplyr)
library(tidyr)
library(JASPAR2020)
library(TFBSTools)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(AcidExperiment)
#library(AcidBase)
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
##Myeloid:
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=unique(combined$celltype2)
for(i in 1:length(cell)){
	print(i)
sub=subset(combined,celltype2==cell[i])
DefaultAssay(sub) <- 'ATAC'
# call peaks for each celltype in myeloid using MACS2
peaks <- CallPeaks(sub, group.by = "AE")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks1 <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks1 <- subsetByOverlaps(x = peaks1, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(sub),
  features = peaks1,
  cells = colnames(sub)
)

fragpath <- Fragments(sub)
# create a new assay using the MACS2 peak set and add it to the Seurat object
sub[["Peaks3"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(sub) <- "Peaks3"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
)
motif.mat <- CreateMotifMatrix(
  features = StringToGRanges(rownames(sub)),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
# add motif information
sub <- AddMotifs(
  object = sub,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
sub <- RegionStats(
  object = sub,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)
#RunChromVAR to caculate the motif activities
sub <- RunChromVAR(
  object = sub,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  motif.matrix = motif.mat,
  new.assay.name = "chromvar"
)

saveRDS(sub,paste0("subset_",cell[i],"_new.rds"))
}

##3 use FigR to get DORC gene
library(FigR)
library(SimBu)
library(EnsDb.Hsapiens.v86)
library(stringr)
library(plyranges)
library(networkD3)
library(htmltools)
library(ggrepel)
require(ggplot2)
require(ggrastr)
require(BuenColors)
library(Seurat) 
library(Signac)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=unique(combined_M$celltype2)
for(i in 1:length(cell)){
print(i)
MM=readRDS(paste0("subset_",cell[i],"_new.rds"))
DefaultAssay(MM) <- "Peaks3"
MM<- RunTFIDF(MM,assays="Peaks3")
counts <- MM@assays$Peaks3@data
counts <- counts[Matrix::rowSums(counts)!=0,]
colData <- DataFrame(MM@meta.data)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
AllPeak <- ClosestFeature(MM, regions = rownames(counts),annotation = annotation,sep = c(':', '-'))
Granges <- data.frame(peak=AllPeak$query_region)
Granges <- str_split_fixed(Granges$peak, ":", 2) 
Granges <- data.frame(chr=Granges[,1],range=Granges[,2])
Granges2 <- str_split_fixed(Granges$range, "-", 2)
Granges <- data.frame(chr=Granges[,1],start=Granges2[,1],end=Granges2[,2])
head(Granges)
Granges <- cbind(AllPeak,Granges)
Granges$start <- as.numeric(Granges$start)
Granges$end <- as.numeric(Granges$end)
peaks_gr <- as_granges(Granges,seqnames = chr)
ATAC.SE <- SummarizedExperiment(assays = list(counts=counts), rowRanges = peaks_gr,colData = colData)

##RNA 
DefaultAssay(MM) <- 'RNA' 
MM <- SCTransform(MM, vst.flavor = "v2",method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE) 
DefaultAssay(MM) <- 'SCT' 
rnaMat= MM@assays$SCT@data
rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
#Peak-gene association by gene's TSS window and peak accessibility(mean-centered) 
cisCor <- runGenePeakcorr(ATAC.se = ATAC.SE,
                           RNAmat = rnaMat,
                           genome = "hg38",
                           nCores = 4, 
                           normalizeATACmat=FALSE,
                           windowPadSize=200000,
                           p.cut=NULL,
                           n_bg = 100)
# Filter peak-gene correlations by p-value                    
cisCor.filt <- cisCor %>% filter(pvalZ <= 0.01)
write.table(cisCor.filt,paste0(cell[i],"_cisCor_filter_200kb.txt"),sep="\t",quote=F)
}


###4 DEG-DORC gene
library(Seurat)
library(Signac)
library(ggplot2)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=unique(combined_M$celltype2)
RES_all=c()
for(i in 1:length(cell)){
  DEG=read.csv(paste0("DEGs_",cell[i],"_res.csv"))
  DORC=read.table(paste0(cell[i],"_cisCor_filter_200kb.txt"),sep="\t",head=T)
  DEG_gene=DEG[which(abs(DEG$log2FoldChange)>0.5 & DEG$padj<0.01),1]
  length(DEG_gene)
  DEG_DORC=intersect(DORC$Gene,DEG_gene)
  if(length(DEG_DORC)>0){
  CIS_RES=c()
  DORC_res=DORC[which(DORC$Gene%in%DEG_DORC),]
  DEG_res=DEG[which(DEG$gene%in%DEG_DORC),]
  colnames(DEG_res)[1]="Gene"
  RES=c()
  for(j in 1:length(DEG_DORC)){
    DEG1=DEG_res[DEG_res$Gene==DEG_DORC[j],c(1,3,7,8)]
    DORC1=DORC_res[DORC_res$Gene==DEG_DORC[j],c(2,4,5)]
    res=cbind(cell[i],DEG1,DORC1)
    RES=rbind(RES,res)
  }
RES_all=rbind(RES_all,RES)
  }
}
colnames(RES_all)=c("Celltype","Gene","log2FoldChange","DEG_padj","threshold","Peak","cis_rObs","cis_Pvalz")
write.table(RES_all,"Mcell_4type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",quote=F,row.names=F,col.names=T)


###Plot DEG gene expression and peak
# define zscore function and apply it to every celltype
zscore_matrix=function(x){
res=c()
for(n in 1:dim(x)[1]){
x_line=x[n,]
x_std <- sd(x_line)*sqrt((length(x_line)-1)/(length(x_line)))
x_mean <- mean(x_line)
x_zscore=(x_line - x_mean) / x_std
res=rbind(res,x_zscore)
}
return(res)
}

library(Rmagic)
library(ggplot2)
library(stringr)
cis=read.table("Mcell_4type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
#order DEG gene by it's counts
DGene=unique(cis$Gene)
R_cc=c()
for(n in 1:length(DGene)){
  cc=length(unique(cis[cis$Gene==DGene[n],1]))
  r_cc=cbind(DGene[n],cc,unique(cis[cis$Gene==DGene[n],5]))
  R_cc=rbind(R_cc,r_cc)
}
R_cc=as.data.frame(R_cc)
R_cc_up=R_cc[which(R_cc[,3]=="UP"),]
R_cc_down=R_cc[which(R_cc[,3]=="Down"),]
R_cc_up_order=R_cc_up[rev(order(as.numeric(R_cc_up[,2]))),]
R_cc_down_order=R_cc_down[rev(order(as.numeric(R_cc_down[,2]))),]
#inter=intersect(R_cc_up[,1],R_cc_down[,1])
#R_cc_up_order[inter]
DGene=unique(c(R_cc_up_order[,1],R_cc_down_order[,1]))#154

#order cell
cell=c("CM","ncM","cDC","pDC")         
res_RNA_all=c()
res_peak_all=c()
col_type_1=c()
for(i in 1:3){
com_sub=readRDS(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/subset_",cell[i],"_new.rds"))
print(cell[i])
##RNA exp,z-score every celltype's gene expression
DefaultAssay(com_sub) = 'RNA'
rnaMat=GetAssayData(com_sub,slot="data")
rnaMat=as.matrix(rnaMat)
rnaMat=rnaMat[DGene,]
rnaMat_zs=zscore_matrix(rnaMat)
rnaMat_Pat_res=c()
rnaMat_Nor_res=c()
peakMat_Pat_res=c()
peakMat_Nor_res=c()
for(n in 1:20){
SS_Pat=sample(1:length(which(com_sub$AE=="Patient")),100)##pDC换成10
SS_Nor=sample(1:length(which(com_sub$AE=="Normal")),100)
rnaMat_Pat=rnaMat_zs[,which(com_sub$AE=="Patient")]
rnaMat_Pat_sample=rnaMat_Pat[,SS_Pat]
rnaMat_Nor=rnaMat_zs[,which(com_sub$AE=="Normal")]
rnaMat_Nor_sample=rnaMat_Nor[,SS_Nor]
rnaMat_Pat_sample_mean=apply(rnaMat_Pat_sample,1,mean)
rnaMat_Nor_sample_mean=apply(rnaMat_Nor_sample,1,mean)
rnaMat_Pat_res=cbind(rnaMat_Pat_res,rnaMat_Pat_sample_mean)
rnaMat_Nor_res=cbind(rnaMat_Nor_res,rnaMat_Nor_sample_mean)
}
rnaMat_sample=cbind(rnaMat_Pat_res,rnaMat_Nor_res)
res_RNA_all=cbind(res_RNA_all,rnaMat_sample)
col_type=cbind(rep(as.character(cell[i]),40),c(rep("Patient",20),rep("Normal",20)))
col_type_1=rbind(col_type_1,col_type)
}
rownames(res_RNA_all)=DGene
write.table(res_RNA_all,"DORC_RNA_exp_zscore_everycelltype_res_sampling.txt",sep="\t",quote=F)
write.table(col_type_1,"zscore_everycelltype_matrix_column_type_sampling.txt",sep="\t",quote=F,row.names=F,col.names=F)

##plot
data_rna=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/DORC_RNA_exp_zscore_everycelltype_res_sampling.txt",sep="\t",head=T)
column_name=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/zscore_everycelltype_matrix_column_type_sampling.txt",sep="\t",head=F)

##order gene
cis=read.table("Mcell_4type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
DEG_gene=rownames(data_rna)
DEG_gene1=DEG_gene[DEG_gene%in%cis[which(cis$Celltype=="CM"&cis$threshold=="UP"),2]]
DEG_gene2=DEG_gene[DEG_gene%in%cis[which(cis$Celltype=="ncM"&cis$threshold=="UP"),2]]
DEG_gene3=DEG_gene[DEG_gene%in%cis[which(cis$Celltype=="cDC"&cis$threshold=="UP"),2]]
DEG_gene4=DEG_gene[DEG_gene%in%cis[which(cis$threshold=="Down"),2]]
DEG_gene12=intersect(DEG_gene1,DEG_gene2)
DEG_gene11=setdiff(DEG_gene1,DEG_gene2)
DEG_gene22=setdiff(DEG_gene2,DEG_gene1)
DEG_gene123=intersect(DEG_gene12,DEG_gene3)
DEG_new=unique(c(DEG_gene123,DEG_gene12,DEG_gene11,DEG_gene22,DEG_gene3,DEG_gene4))

data_rna1=data_rna[DEG_new,]
colnames(column_name)=c("Celltype","Group")
data_rna1[is.na(data_rna1)]=0
cell_order=c(which(column_name$Group=="Patient"),which(column_name$Group=="Normal"))
data_rna1=data_rna1[,cell_order]
column_name=column_name[cell_order,]
#rownames(column_name)=1:length(column_name[,1])
 # creat colours for each group 
group_color=paletteer::paletteer_d("lisa::OskarSchlemmer")
cell_color=paletteer::paletteer_d("ggpomological::pomological_palette")
annoColor <- list(Group=c("Normal"=group_color[1], "Patient"=group_color[5]),
                  Celltype=c("CM"=cell_color[1],"ncM"=cell_color[2],"cDC"=cell_color[3]))
#data_rna1=t(data_rna)
colnames(data_rna1)=rownames(column_name)
#rownames(data32)=rownames(ann_row)
library(pheatmap)
pdf("Fig5B_DORC_DEG_RNA_pheatmap.pdf",height =8,width=5)#Fig6C
pheatmap(data_rna1,    #labels_row=DGene, 
color = colorRampPalette(paletteer::paletteer_d("khroma::BuRd"))(100),
breaks=(c(seq(-0.5, 0, length.out = 50),seq(0.001, 0.5, length.out = 50))),
cluster_rows = FALSE,   cluster_cols = FALSE,
show_rownames = TRUE,show_colnames = FALSE, #annotation_row = ann_row,
annotation_col = column_name,  annotation_colors = annoColor,
border_color = NA,      fontsize = 6,      scale = "none",#file="DORC_DEG_RNA_pheatmap_order.pdf"
) 
dev.off()

2.Fig 5B
#####将三个celltype上调的的基因富集的功能拿出来做交集，展示一些top的功能在C图里面
#setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5/")
#对所有上调的108个基因做富集分析,展示其中部分term
##Metascape 做富集分析
#Supplementary Table2




3.Fig 5C 受体蛋白的表达情况
#因为Myelid作为固有免疫，是免疫反应的初始应答过程，固有免疫细胞识别抗原进行吞噬或者识别，能够吞噬抗原或者受到抗原激活后激活其他免疫细胞
###这个过程中涉及的模式识别受体（胞膜或者胞质）包括几种类别
##分别是Toll样受体(TLR)，NOD样受体(NLR)，RIG样受体（RLR)和胞质DNA感知蛋白（ cGAS）
##TLR
TLR=c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TLR10")
NOD=c("NLRP1","NLRP3","NLRC4","NOD1","NOD2","CIITA","NAIP",#NLRA--CIITA,"NLRB"--NAIP
      "NLRA1","NLRA2","NLRA3","NLRA4","NLRB1","NLRC1","NLRC2","NLRC3","NLRC5","NLRD1","NLRF1",
        "NLRF2","NLRF3","NLRF4","NLRF5","NLRF6","NLRG1","NLRG2","NLRH1","NLRH2","NLRH3","NLRH4",
        "NLRX1","NLRY1","NLRY2","NLRX3","NLRX4","NLRX5","NLRX6",
        "NLRP2","NLRP5","NLRP6","NLRP7","NLRP8","NLRP9","NLRP10","NLRP11","NLRP12",
        "NLRP13","NLRP14","NLRP15","NLRP16","NLRP17","NLRP18","NLRP19","NLRP20")
RLR_cGAS=c("DDX58","IFIH1","DDX60L",
           "CGAS","TMEM173","MPYS","ZBP1","AIM2")#"STING"--TMEM173
           #DNA依赖型干扰素调节因子激活蛋白(DAI)--ZBP1
           #黑色素瘤缺失因子--AIM2

##检验基因名
DefaultAssay(combined_M)="RNA"
M_RNA=GetAssayData(combined_M)
TLR_1=TLR[TLR%in%rownames(M_RNA)]
NOD_1=NOD[NOD%in%rownames(M_RNA)]
# [1] "NLRP1"  "NLRP3"  "NLRC4"  "NOD1"   "NOD2"   "CIITA"  "NAIP"   "NLRC3" 
# [9] "NLRC5"  "NLRX1"  "NLRP2"  "NLRP5"  "NLRP6"  "NLRP7"  "NLRP8"  "NLRP9" 
# [17] "NLRP10" "NLRP11" "NLRP12" "NLRP13" "NLRP14"
RLR_cGAS_1=RLR_cGAS[RLR_cGAS%in%rownames(M_RNA)]
gene1=c(TLR,NOD_1,RLR_cGAS_1)

#test Patient vs Normal pvalue
combined_M$cellAE=paste0(combined_M$celltype2,"_",combined_M$AE)
Idents(combined_M)="cellAE"
DefaultAssay(combined_M)="RNA"
RNA_data=GetAssayData(combined_M,slot="data")
RNA_data1=RNA_data[gene1,]
cell=unique(combined_M$celltype2)
RES=c()
for(i in 1:length(gene1)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[gene1[i],which(combined_M$celltype2==cell[n]&combined_M$AE=="Patient")]
  exp_Nor=RNA_data1[gene1[i],which(combined_M$celltype2==cell[n]&combined_M$AE=="Normal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  res=cbind(gene1[i],as.character(cell[n]),Test_p,FC)
  RES=rbind(RES,res)
}
}
colnames(RES)=c("Gene","Cell","Pvalue","FC")
write.table(RES,"Fig5C_Receptor_gene_Pat_vs_Nor_test.txt",sep="\t",quote=F,row.names=F)


gene=c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8",
"TLR9","TLR10","NLRP1","NLRP3","NLRC4","NOD1","NOD2","CIITA","NAIP","NLRC3","NLRC5","NLRX1","NLRP2",
"NLRP5","NLRP6","NLRP7","NLRP8","NLRP9","NLRP10","NLRP11","NLRP12","NLRP13","NLRP14","DDX58","IFIH1","DDX60L",
"CGAS","MPYS","ZBP1","AIM2")
RES=read.table("Fig5C_Receptor_gene_Pat_vs_Nor_test.txt",sep="\t",head=T)
DefaultAssay(combined_M)="RNA"
M_RNA=GetAssayData(combined_M)
for(i in 1:length(RES[,1])){
  RES[i,5]=mean(M_RNA[RES[i,1],which(combined_M$celltype2==RES[i,2])])
  RES[i,6]=length(which(M_RNA[RES[i,1],which(combined_M$celltype2==RES[i,2])]>0))/length(M_RNA[RES[i,1],which(combined_M$celltype2==RES[i,2])])
}
colnames(RES)[c(5,6)]=c("exp_mean","exp_Ratio")
write.table(RES,"Fig5C_Receptor_gene_Pat_vs_Nor_test.txt",sep="\t",quote=F,row.names=F)

RES=read.table("Fig5C_Receptor_gene_Pat_vs_Nor_test.txt",sep="\t",head=T)
RES=RES[which(RES$Gene%in%c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8",
"TLR9","TLR10","NLRP1","NLRP2","NLRP3","NLRP5","NLRP6","NLRP7","NLRP8","NLRP9","NLRP10","NLRP11","NLRP12",
"NOD1","NOD2","NLRC3","NLRC4","NLRC5","NLRX1","CIITA","NAIP",
"DDX58","IFIH1","DDX60L",
"CGAS")),]
pdf("Fig5C_Receptor_gene_Pat_vs_Nor_test_new.pdf",width=5,height=8)
ggplot(RES,aes(x=Cell,y=factor(Gene,levels=c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8",
"TLR9","TLR10","NLRP1","NLRP2","NLRP3","NLRP5","NLRP6","NLRP7","NLRP8","NLRP9","NLRP10","NLRP11","NLRP12",
"NOD1","NOD2","NLRC3","NLRC4","NLRC5","NLRX1","CIITA","NAIP",
"DDX58","IFIH1","DDX60L",
"CGAS")),
  size=exp_Ratio,color=log2(FC)))+geom_point(alpha=1) + 
	labs(x= "TF activity(log2FC)") +ylab("TF expression(Ratio)") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
  theme_bw()+
  scale_size(breaks = c(0,0.1,0.2,0.3),range = c(-1,10))
dev.off()

#key_gene=unique(RES[which(RES$exp_Ratio>0.1&RES$exp_mean>0.2&RES$FC>1&RES$Pvalue<0.01),1])
#[1] "TLR2"   "NLRP1"  "NLRP3"  "NAIP"   "NLRC5"  "NLRP12"
#key_gene=unique(RES[which(RES$exp_Ratio>0.2&RES$exp_mean>0.2&RES$FC>1&RES$Pvalue<0.01),1])


4.Fig 5D vlnplot of same receptor genes
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
DefaultAssay(combined_M)="RNA"
Idents(combined_M)="celltype2"
cell=c("CM","ncM","cDC","pDC")
names(cell) <- levels(combined_M)
combined_M <- RenameIdents(combined_M, cell)
test_sign=list(c("Patient","Normal"))
gene=c("TLR2","NLRP1","NLRP3","NLRP12","NAIP","NLRC5")
pdf("Fig5D_TLR2_ex_gene_exp_mean_vlnplot_RNA.pdf",width=6,height=4)
for(m in 1:length(gene)){
p=VlnPlot(combined_M, features = gene[m], 
        split.by = 'AE', pt.size = 0.0, 
        combine = T, split.plot = TRUE,##Seurat<'4.3.0.1',可以是4.1.1
        cols =c("#437DBFFF","#AA6F9EFF")) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
               position = position_dodge(width = .70)) +   
       xlab("Cell") + ylab("Gene Expression") + 
       #stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       #stat_compare_means( method = "wilcox.test") +
       theme(text = element_text(size = 8, angle = 0))+
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
      print(p)
}
dev.off()



5. Fig 5E
###TNF pathway gene expression
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(ggplot2)
library(devtools)
library(ggpubr)
library(dplyr)
TNF=c("TRADD","TRAF2",
       "MAP3K5","MAP3K4","MAP3K7","MAPK8","MAPK9","JUN",#AP-1
      "RIPK1","MAP3K3","MAP3K6", "MAPK1","MAPK3",#MAPK
      "MAP3K14","IKBKB","CHUK","IKBKB","RELA","REL","NFKB1","NFKB2"#NFKB
      )
      #IKKa-CHUK,IKKb-IKBKB,RIP-RIPK1,MAPK-MAPK1/MAPK3;
      #ASK1-MAP3K5;JNK-MAPK8/MAPK9/MAPK10,AP1-JUN
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
DefaultAssay(combined)="RNA"
RNA_data=GetAssayData(combined,slot="data")
RNA_data1=RNA_data[TNF,]
#TNF%in%rownames(RNA_data)
#unique(combined$celltype2)
cell=c("B naive","B memory","B plasma",
"CD4 naive","CD4 Tm", "CD4 Treg" ,"CD8 naive","CD8 Tem","CD8 Tcm","CD8 MAIT","CD8 Treg" ,"NK")       
RES=c()
for(i in 1:length(TNF)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[TNF[i],which(combined$celltype2==cell[n]&combined$AE=="Patient")]
  exp_Nor=RNA_data1[TNF[i],which(combined$celltype2==cell[n]&combined$AE=="Normal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  Ratio_Pat=length(which(exp_Pat>0))/length(exp_Pat)
  mean_Pat=mean(exp_Pat)
  Ratio_Nor=length(which(exp_Nor>0))/length(exp_Nor)
  mean_Nor=mean(exp_Nor)

  res1=cbind(TNF[i],cell[n],Test_p,FC,Ratio_Pat,mean_Pat,"Patient")
  res2=cbind(TNF[i],cell[n],Test_p,FC,Ratio_Nor,mean_Nor,"Normal")
  RES=rbind(RES,res1,res2)
}
}
colnames(RES)=c("Gene","Cell","Pvalue","FC","Ratio_exp","Mean_exp","Type")
write.table(RES,"Fig5E_TNF_pathway_gene_Pat_vs_Nor_test.txt",sep="\t",quote=F,row.names=F)


RES=read.table("Fig5E_TNF_pathway_gene_Pat_vs_Nor_test.txt",sep="\t",head=T)##去掉NA,inf改成mean(Patient)
RES$CELL=paste0(RES$Cell,"-",RES$Type)
pdf("Fig5E_TNF_pathway_gene_Pat_vs_Nor_test_allCell.pdf",width=12,height=5)
ggplot(RES,aes(x=factor(CELL,levels=c("CM-Patient","CM-Normal","ncM-Patient","ncM-Normal",
"cDC-Patient","cDC-Normal","pDC-Patient","pDC-Normal",
"B naive-Patient","B naive-Normal","B memory-Patient","B memory-Normal","B plasma-Patient","B plasma-Normal",
"CD4 naive-Patient","CD4 naive-Normal","CD4 Tm-Patient","CD4 Tm-Normal","CD4 Treg-Patient","CD4 Treg-Normal",
"CD8 naive-Patient","CD8 naive-Normal","CD8 Tcm-Patient","CD8 Tcm-Normal","CD8 Tem-Patient","CD8 Tem-Normal",
"NK-Patient","NK-Normal","CD8 MAIT-Patient","CD8 MAIT-Normal","CD8 Treg-Patient","CD8 Treg-Normal")),
  y=factor(Gene,levels=c("TNFSF13B","TNFSF10","TNFRSF13B","TNFRSF13C","TNFRSF17","TNFRSF10A")),
  size=Ratio_exp,color=log2(FC)))+geom_point()+
	labs(x= "cell") +ylab("Gene") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
  theme_bw()+
  scale_size(breaks = c(0,0.05,0.1,0.15,0.2),range = c(0,10))+
  theme(axis.text.x = element_text(size = 10, angle = 90),
axis.text.y = element_text(size = 6, angle = 0))
dev.off()











































###先看一下分泌蛋白（细胞因子）的表达情况
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/")
IL=c("IL1A","IL1B","IL2","IL3","IL4","IL5","IL6","IL7","IL8","IL9","IL10","IL11","IL12A",
    "IL12B","IL13","IL15","IL16","IL17A","IL17B","IL17C","IL17D","IL17E","IL17F","IL18","IL19",
    "IL20","IL21","IL22","IL23A","IL24","IL25","IL26","IL27","IL28A","IL28B","IL29","IL31",
    "IL32","IL33","IL34","IL35","IL36A","IL36B","IL36G","IL37")
TNF=c("TNF","TNFAIP3","TNFRSF1A","TNFRSF1B","TNFRSF4","TNFRSF8","TNFRSF9","TNFRSF10A",
     "TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF11A","TNFRSF11B","TNFRSF12A","TNFRSF12B",
     "TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF19","TNFRSF21","TNFRSF22",
     "TNFRSF25","TNFRSF6B","TNFSF4","TNFSF5","TNFSF8","TNFSF9","TNFSF10","TNFSF11","TNFSF12",
     "TNFSF13B","TNFSF13","TNFSF14","TNFSF15","TNFSF18","TNFSF13B","TNFSF13","TNFSF14",
     "TNFSF15","TNFSF18","TNFSF20","TNFSF25","TNFSF11")
IFN=c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNB1","IFNG","IFNE","IFNK","IFNW1","IFNW2")
GF=c("EGF","FGF1","FGF2","FGF3","FGF4","FGF5","FGF6","FGF7","FGF8","FGF9","FGF10","FGF11",
     "FGF12","FGF13","FGF14","FGF16","FGF17","FGF18","FGF19","FGF20","FGF21","FGF22","FGF23",
     "PDGF","TGFB1","TGFB2","TGFB3","HGF","VEGFA","VEGFB","VEGFC","VEGFD","NGF")
CCL=c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL6","CCL7","CCL8","CCL9","CCL10","CCL11",
     "CCL12","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21",
     "CCL22","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CXCL1","CXCL2","CXCL3",
     "CXCL4","CXCL5","CXCL6","CXCL7","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
DefaultAssay(combined)="RNA"
IL=unique(IL)
#dotplot的改造  
Idents(combined)="celltype2"
unique(combined$celltype2)
pdf("Fig6_Dotplot_IL_gene_RNA.pdf",width=7,height=8)
DotPlot(combined, features = IL,dot.scale = 8)+
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
)#"Progenitor",
dev.off()

IL=c("IL1B","IL2","IL2A","IL6","IL7","IL10","IL12RB1","IL13","IL15","IL16",
     "IL17D","IL18","IL19","IL23A","IL24","IL26","IL32")
pdf("Fig6_Dotplot_IL_gene_RNA_new.pdf",width=6,height=6)
DotPlot(combined, features = IL,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,3,5,10,20, 30),range = c(-1,8))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)#"Progenitor",
dev.off()


TNF=unique(TNF)
pdf("Fig6_Dotplot_TNF_gene_RNA.pdf",width=7,height=8)
DotPlot(combined, features = TNF,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,3,5,10,20,30,40),range = c(-1,8))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)
dev.off()

TNF=c("TNF","TNFAIP3","TNFRSF1A","TNFRSF1B","TNFRSF8","TNFRSF10A",
     "TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF11A",
     "TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF21",
     "TNFRSF25","TNFSF8","TNFSF10",
     "TNFSF13B","TNFSF13","TNFSF14","TNFSF15")
pdf("Fig6_Dotplot_TNF_gene_RNA_new.pdf",width=7,height=8)
DotPlot(combined, features = TNF,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,3,5,10,20,30,40),range = c(-1,8))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)
dev.off()

IFN=unique(IFN)
pdf("Fig6_Dotplot_IFN_gene_RNA.pdf",width=7,height=8)
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
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)#"Progenitor",
dev.off()


GF=unique(GF)
pdf("Fig6_Dotplot_GF_gene_RNA.pdf",width=7,height=8)
DotPlot(combined, features = GF,dot.scale = 8)+
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
)#"Progenitor",
dev.off()

CCL=unique(CCL)
pdf("Fig6_Dotplot_CCL_gene_RNA.pdf",width=7,height=8)
DotPlot(combined, features = CCL,dot.scale = 8)+
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
)#"Progenitor",
dev.off()

IFN_GF_CCL=c("IFNG","TGFB1","TGFB2","HGF","VEGFA","VEGFB","CCL3","CCL4","CCL5","CXCL8","CCL28","CXCL13")
pdf("Fig6_Dotplot_IFN_GF_CCL_gene_RNA_new.pdf",width=6,height=5)
DotPlot(combined, features = IFN_GF_CCL,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,5,10, 20, 30,40),range = c(-1,10))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)#"Progenitor",
dev.off()



IL=c("IL1B","IL2","IL2A","IL6","IL7","IL10","IL12RB1","IL13","IL15","IL16",
     "IL17D","IL18","IL19","IL23A","IL24","IL26","IL32")
IFN_GF_CCL=c("IFNG","TGFB1","TGFB2","HGF","VEGFA","VEGFB","CCL3","CCL4","CCL5","CXCL8","CCL28","CXCL13")
TNF=c("TNF","TNFAIP3","TNFRSF1A","TNFRSF1B","TNFRSF8","TNFRSF10A",
     "TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF11A",
     "TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF21",
     "TNFRSF25","TNFSF8","TNFSF10",
     "TNFSF13B","TNFSF13","TNFSF14","TNFSF15")
other=c("NLRP3","CD40L","CD40","CD40LG")
gene_new=c(IL,IFN_GF_CCL,TNF,other)
pdf("Fig6_Dotplot_all_immune_gene_RNA_new.pdf",width=8,height=10)
DotPlot(combined, features = gene_new,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,5,10, 20, 30,40),range = c(-1,10))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive","B memory","B plasma","CD4 naive","CD4 Tm",
  "CD4 Treg","CD8 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","NK","CM",
  "ncM","cDC","pDC")
)#"Progenitor",
dev.off()

###在Myeloid中表达的因子单独看差异性

M_gene=c("IL1B","IL15","IL18","HGF","CXCL8","TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF21",
     "TNFSF10","TNFSF13B","TNFSF13","TNFSF15","NLRP3")
com_M=subset(combined,celltype2%in%c("cDC","CM","ncM","pDC"))
com_M$Cell_AE=paste0(com_M$celltype2,"_",com_M$AE)
Idents(com_M)="Cell_AE"
pdf("Fig6_Myeloid_Dotplot_all_immune_gene_RNA_new.pdf",width=6,height=6)
DotPlot(com_M, features = M_gene,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,5,10, 20, 30,40),range = c(-1,10))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("CM_Normal","CM_Patient","ncM_Normal","ncM_Patient",
                            "cDC_Normal","cDC_Patient","pDC_Normal","pDC_Patient")
)
dev.off()

B_gene=c("CD40","IL7","IL6","TNFRSF17","TNFRSF13C","TNFRSF13B","TNFRSF10A")
com_B=subset(combined,celltype2%in%c("B memory","B naive","B plasma"))
com_B$Cell_AE=paste0(com_B$celltype2,"_",com_B$AE)
Idents(com_B)="Cell_AE"
pdf("Fig6_BCell_Dotplot_all_immune_gene_RNA_new.pdf",width=6,height=6)
DotPlot(com_B, features = B_gene,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,5,10, 20, 30,40),range = c(-1,10))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("B naive_Normal","B naive_Patient","B memory_Normal","B memory_Patient",
                            "B plasma_Normal","B plasma_Patient")
)
dev.off()

T_gene=c("CD40LG","CCL4","CCL5","TNFSF14","TNFSF8","TNFRSF25","TNFRSF10A","TGFB1","IFNG","IL32")
com_T=subset(combined,celltype2%in%c("CD4 naive","CD4 Tm","CD4 Treg","CD8 MAIT","CD8 naive","CD8 Tcm",
                                     "CD8 Tem","CD8 Treg","NK"))
com_T$Cell_AE=paste0(com_T$celltype2,"_",com_T$AE)
Idents(com_T)="Cell_AE"
pdf("Fig6_TCell_Dotplot_all_immune_gene_RNA_new.pdf",width=9,height=6)
DotPlot(com_T, features = T_gene,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,5,10, 20, 30,40),range = c(-1,10))+
labs(y="Cell Type",x="Marker Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits = c("CD4 naive_Normal","CD4 naive_Patient","CD4 Tm_Normal","CD4 Tm_Patient",
                            "CD4 Treg_Normal","CD4 Treg_Patient","CD8 MAIT_Normal","CD8 MAIT_Patient",
                            "CD8 naive_Normal","CD8 naive_Patient","CD8 Tcm_Normal","CD8 Tcm_Patient",
                            "CD8 Tem_Normal","CD8 Tem_Patient","CD8 Treg_Normal","CD8 Treg_Patient",
                            "NK_Normal","NK_Patient")
)
dev.off()


#####某些关键TNF因子的track
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
M_gene=c("IL1B","IL15","IL18","HGF","CXCL8","TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFRSF21",
     "TNFSF10","TNFSF13B","TNFSF13","TNFSF15","NLRP3")
com_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")   
DefaultAssay(com_M)="Peaks1"
com_M <- RegionStats(com_M, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(com_M) <- 'SCT' 
SCT=GetAssayData(com_M,slot="data")
Gene=M_gene[which(M_gene%in%rownames(SCT))]
com_M <- LinkPeaks(
  object = com_M,
  peak.assay = "Peaks1",
  expression.assay = "SCT",
  genes.use = Gene,
  pvalue_cutoff = 0.01,
  score_cutoff = 0
  )
com_M$cell_Label=paste0(com_M$celltype2,"_",com_M$AE)
Idents(com_M) <- "cell_Label"
levels(com_M) <- c("CM_Normal","CM_Patient","ncM_Normal","ncM_Patient",
                            "cDC_Normal","cDC_Patient","pDC_Normal","pDC_Patient"
                   )
DefaultAssay(com_M)="Peaks1"

pdf("Fig6_Myeloid_immune_gene_track_200kb.pdf",height=6,width=10)
for(i in c(1:length(Gene))){
p=CoveragePlot(
   object = com_M,
   region = Gene[i], 
   #region.highlight =ranges.show,
   #ranges = Peaks1,
   features = Gene[i],
   expression.assay = "SCT",
   extend.upstream = 200000,
   extend.downstream = 200000,
   #ymax = 20,##change normalized signal range
   links=T,
   window = 200
 )
 print(p)
}
 dev.off()

com_ncM=readRDS("subset_ncM_new.rds")
DefaultAssay(com_ncM)="Peaks3"
com_ncM <- RegionStats(com_ncM, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(com_ncM) <- 'SCT' 
SCT=GetAssayData(com_ncM,slot="data")
Gene=M_gene[which(M_gene%in%rownames(SCT))]
com_ncM <- LinkPeaks(
  object = com_ncM,
  peak.assay = "Peaks3",
  expression.assay = "SCT",
  genes.use = Gene,
  pvalue_cutoff = 0.01,
  score_cutoff = 0
  )
  DefaultAssay(com_ncM)="Peaks3"
  Idents(com_ncM)="AE"
pdf("Fig6_Myeloid_immune_gene_track_200kb_nCM.pdf",height=6,width=10)
for(i in 1:length(Gene)){
p=CoveragePlot(
   object = com_ncM,
   region = Gene[i], 
   #region.highlight =ranges.show,
   #ranges = Peaks1,
   features = Gene[i],
   expression.assay = "SCT",
   extend.upstream = 100000,
   extend.downstream = 100000,
   #ymax = 20,##change normalized signal range
   links=T,
   window = 200
 )
 print(p)
 }
 dev.off()

















###################################################################################################
####subset Monocytes
library(Seurat)
library(ggplot2)
library(dittoSeq)
library(MySeuratWrappers)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig6/Monocyte_subtype/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=table(combined_M$celltype2)
com=subset(combined_M,celltype2=="CM")
DefaultAssay(com) <- "RNA"
com <- FindMultiModalNeighbors(com, reduction.list = list("pca", "integratedLSI"), dims.list = list(1:30, 2:50))
com <- RunUMAP(com, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
com <- FindClusters(com, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(com$seurat_clusters)
table(com@meta.data[,c(48,58)])/apply(table(com@meta.data[,c(48,58)]),1,sum)

pdf("Dimplot_M_clusters.pdf",width=5,height=5)
DimPlot(com, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black")
dev.off()
pdf("Dimplot_M_clusters1.pdf",width=5,height=5)
DimPlot(com, reduction = "wnn.umap", group.by = 'Label', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black")
dev.off()

pdf("Dimplot_M_clusters_1.pdf",width=25,height=5)
DimPlot(com, reduction = "wnn.umap", group.by = 'seurat_clusters', pt.size = 0.2,label = T,
  label.size = 4, label.color = "black",split.by="Label")
dev.off()
####看起来并没有特别大的patient特异性
