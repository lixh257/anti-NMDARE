Fig4
#####
1.Fig 4A
library(Seurat)
library(ggplot2)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4")
combined_T=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")

pdf("Fig4A_Dimplot_T_celltype.pdf",height=5,width=6)
DimPlot(combined_T, reduction = "wnn.umap", group.by = 'celltype2', pt.size = 0.1,label = T,
  label.size = 4, label.color = "black")
dev.off()


2. Fig 4B
##barplot of every individual cell's propotion
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
colnames(combined_T@meta.data)
table(combined_T@meta.data[,c(48,60)])
data=table(combined_T@meta.data[,c(48,60)])/apply(table(combined_T@meta.data[,c(48,60)]),1,sum)
data=as.data.frame(data)

pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#A3D8F4","#8CC0DE","#A1CCD1","#64C9CF","#01A9B4","#176B87")
pdf("Fig4B_barplot_cellratio_of_Tcell_Sample.pdf",width=12,height=8)
ggplot(data,
aes(x=factor(celltype2,levels=c("NK","CD8 Tem","CD8 Tcm","CD8 MAIT","CD8 Treg","CD8 naive","CD4 naive","CD4 Tm","CD4 Treg")),
y=Freq,
fill=factor(Label,levels=c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Patient1","Patient3","Patient4"))))+
geom_bar(stat="identity", width=.8, position = "dodge")+
scale_fill_manual(values=c(pal_Nor,pal_pat),name="Sample")+
theme_bw()+
labs(x="Celltype",y="The proportion of each individual")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()

####Fig 4B 右上角
sample_table <- as.data.frame(table(combined_T@meta.data$AE,combined_T@meta.data$celltype2))
names(sample_table) <- c("Samples","celltype","CellNumber")
sample_table$Samples <- factor(sample_table$Samples, level=c("Normal","Patient"))
sample_table$celltype <- factor(sample_table$celltype, level=c("NK","CD8 Tem","CD8 Tcm","CD8 MAIT","CD8 Treg","CD8 naive","CD4 naive","CD4 Tm","CD4 Treg"))

pal=c("#2CA02CFF" ,"#98DF8AFF","#DBDB8DFF","#E377C2FF","#AEC7E8FF","#F7B6D2FF","#9467BDFF","#7F7F7FFF","#176B87")
pdf("Fig4B_Percentage-celltype-in-sample_group.pdf")                     
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


3-6.Fig4C-F
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
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
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
DefaultAssay(combined) <- 'RNA' 
Idents(combined) <- 'AE' 
cell=unique(combined$celltype2)
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
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
##Tcell:
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4")
comNKT=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
#names(NKT)
cell=unique(comNKT$celltype2)
#[1] CD4 Tm    CD4 naive CD8 naive CD8 Tcm   CD8 MAIT  CD8 Tem   NK Treg
for(i in 1:length(cell)){
	print(i)
sub=subset(comNKT,celltype2==cell[i])
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
library(Signac)
library(Seurat)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
comNKT=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
#names(NKT)
cell=unique(comNKT$celltype2)
for(i in 1:9){
print(i)
NKT=readRDS(paste0("subset_",cell[i],"_new.rds"))
DefaultAssay(NKT) <- "Peaks3"
NKT<- RunTFIDF(NKT,assays="Peaks3")
counts <- NKT@assays$Peaks3@data
counts <- counts[Matrix::rowSums(counts)!=0,]
colData <- DataFrame(NKT@meta.data)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
AllPeak <- ClosestFeature(NKT, regions = rownames(counts),annotation = annotation,sep = c(':', '-'))
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
DefaultAssay(NKT) <- 'RNA' 
NKT <- SCTransform(NKT, vst.flavor = "v2",method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE,return.only.var.genes = FALSE) 
DefaultAssay(NKT) <- 'SCT' 
rnaMat= NKT@assays$SCT@data
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
library(stringr)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4")
combined_T=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
cell=unique(combined_T$celltype2)
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
write.table(RES_all,"Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",quote=F,row.names=F,col.names=T)



3. Fig 4C
#####统计差异基因的数量并展示其在每种细胞类型中的分布
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
cell=c("CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK") 
res5=c()
for(i in 1:length(cell)){  
#DEG=cis[which(cis$Celltype==cell[i]),]
DEG=read.csv(paste0("DEGs_",cell[i],"_res.csv"))
#order cell
res1=length(unique(DEG[which(DEG$threshold=="UP"),1]))
res2=length(unique(DEG[which(DEG$threshold=="Down"),1]))
res3=c(cell[i],"Up",res1)
res4=c(cell[i],"Down",res2)
res5=rbind(res5,res3,res4)
}
colnames(res5)=c("Celltype","Regu","number")
write.table(res5,"Fig4C_DEG_number_res.txt",sep="\t",quote=F,row.names=F)

res=read.table("Fig4C_DEG_number_res.txt",sep="\t",head=T)
library(ggplot2)
library(plyr)
library(gridExtra)

pal <- paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20, 3,13,19, 7,15,17,5,14,11, 12,16,18,8)]
p_up=ggplot(subset(res,Regu=="Up"),aes(x=factor(Celltype,levels=c("CD4 naive","CD4 Tm","CD4 Treg",
"CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))

p_down=ggplot(subset(res,Regu=="Down"),aes(x=factor(Celltype,levels=c(
"CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem",
"CD8 MAIT","CD8 Treg","NK")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
## Plutting it together
pdf("Fig4C_DEG_number.pdf",width=12,height=5)
grid.arrange(p_up, p_down, ncol=2)
dev.off()

##cis_peak的数量
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
cell=c("CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK") 
res5=c()
for(i in 1:length(cell)){  
DEP=cis[which(cis$Celltype==cell[i]),]
res1=length(unique(DEP[which(DEP$threshold=="UP"),6]))
res2=length(unique(DEP[which(DEP$threshold=="Down"),6]))
res3=c(cell[i],"Up",res1)
res4=c(cell[i],"Down",res2)
res5=rbind(res5,res3,res4)
}
colnames(res5)=c("Celltype","Regu","number")
write.table(res5,"Fig4C_DEP_number_res.txt",sep="\t",quote=F,row.names=F)

res=read.table("Fig4C_DEP_number_res.txt",sep="\t",head=T)
library(ggplot2)
library(plyr)
library(gridExtra)

pal <- paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20, 3,13,19, 7,15,17,5,14,11, 12,16,18,8)]
p_up=ggplot(subset(res,Regu=="Up"),aes(x=factor(Celltype,levels=c("CD4 naive","CD4 Tm","CD4 Treg",
"CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Peak number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))

p_down=ggplot(subset(res,Regu=="Down"),aes(x=factor(Celltype,levels=c(
"CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem",
"CD8 MAIT","CD8 Treg","NK")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Peak number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
## Plutting it together
pdf("Fig4C_DEP_number.pdf",width=12,height=5)
grid.arrange(p_up, p_down, ncol=2)
dev.off()



4.Fig 4D  get peak enrichment motif
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
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4")
combined_T=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
cell=unique(cis$Celltype)
RES_all=c()
for(i in 1:length(cell)){
  sub=readRDS(paste0("subset_",cell[i],"_new.rds"))
  markerPeaks=cis[which(cis$Celltype==cell[i]&cis$threshold=="UP"),6]
  markerPeaks=gsub(":","-",markerPeaks)
###关注上调的gene-link-peak list
## make GRanger Format of MarkerPeak list
T_Granges <- str_split_fixed(markerPeaks, "-", 2) 
T_Granges <- data.frame(chr=T_Granges[,1],range=T_Granges[,2])
T_Granges2 <- str_split_fixed(T_Granges$range, "-", 2)
T_Granges <- data.frame(chr=T_Granges[,1],start=T_Granges2[,1],end=T_Granges2[,2])
T_markerPeak_site <- cbind(markerPeaks,T_Granges)
T_markerPeak_site$start <- as.numeric(T_markerPeak_site$start)
T_markerPeak_site$end <- as.numeric(T_markerPeak_site$end)
T_markerPeak_site <- makeGRangesFromDataFrame(T_markerPeak_site)
markerpeak_gene <- annotatePeak(T_markerPeak_site, TxDb=txdb,tssRegion=c(-3000, 3000), 
                                addFlankGeneInfo=TRUE, flankDistance=200000, annoDb="org.Hs.eg.db")
#markerpeak_gene
write.csv(data.frame(markerpeak_gene),paste0(cell[i],"_DEG_T-up-cisPeak.csv"),quote=F)

pdf(paste0(cell[i],"_DEG_T-up-markerPeak-distribution.pdf"))
plotAnnoBar(markerpeak_gene)
plotDistToTSS(markerpeak_gene,title="Distribution of transcription MarkerPeak loci \n relative to TSS")
dev.off()

library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(dplyr)
DefaultAssay(sub)="Peaks3"
sub <- RegionStats(sub, genome = BSgenome.Hsapiens.UCSC.hg38)
T_markerPeak_motif <- FindMotifs(sub, features = markerPeaks)
# Rank motifs by p-value and fold enrichment
T_markerPeak_motif_rank <- T_markerPeak_motif[order(T_markerPeak_motif$pvalue, T_markerPeak_motif$fold.enrichment), ]
write.table(T_markerPeak_motif_rank,"1-T_markerPeak_motif_rank.txt",sep="\t",quote=F)
#extract motif name
#Motif_TF_name <- T_markerPeak_motif_rank[,c(1,8)]
#write.table(Motif_TF_name,"Motif_TF_name.txt",sep="\t",quote=F)

##### Specific motif by differential chromvar score #####
DefaultAssay(sub) <- 'chromvar'
Idents(sub) <- "AE"
T_chromvarMotif_vsNormal <- FindMarkers(
  object = sub,
  ident.1="Patient",
  ident.2="Normal",
  min.pct = 0.1,
  logfc.threshold=0,  # Default 0.25
  #only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
T_chromvarMotif_vsNormal$motif=rownames(T_chromvarMotif_vsNormal)
Motif_res1=merge(T_markerPeak_motif_rank,T_chromvarMotif_vsNormal,by = "motif")
Motif_res1=Motif_res1[order(Motif_res1$pvalue),]
Motif_res1=Motif_res1[,c(8,1,4,6,9,14,11)]
colnames(Motif_res1)=c("motif.name","motif","enrich_percent.observed",
                       "enrich_fold.enrichment","enrich_padj","chromvar_padj","chromvar_avg_diff")

#add expression of the motif
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "AE"
aa=GetAssayData(sub)
aa_mean_P=apply(aa[,which(sub$AE=="Patient")],1,mean)
aa_mean_N=apply(aa[,which(sub$AE=="Normal")],1,mean)
exp_mean=rbind(aa_mean_P,aa_mean_N)
exp_mean=t(exp_mean)
rname=rownames(exp_mean)
exp_mean=as.data.frame(exp_mean)
colnames(exp_mean)=c("Patient","Normal")
rownames(exp_mean)=rname

###motif.name
Motif=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/Motif_TF_name.txt",sep="\t",head=T)###手动修改这个cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5_Bcell/ /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2
#setdiff(motif$motif.name.1,rownames(exp_mean))
Motif_res1$RNA_log2FC=rep(0,length(Motif_res1[,1]))
Motif_res1$Pat_mean_RNA=rep(0,length(Motif_res1[,1]))
Motif_res1$Pat_RNA_ratio=rep(0,length(Motif_res1[,1]))
for(m in 1:length(Motif_res1[,1])){
  tf_gene=Motif[which(Motif[,2]==Motif_res1[m,2]),4]
  RNA_FC=(exp_mean[tf_gene,1]+0.001)/(exp_mean[tf_gene,2]+0.001)
  RNA_FC_m=sum(RNA_FC)/length(RNA_FC)
  Motif_res1$RNA_log2FC[m]=log2(RNA_FC_m)

  RNA_mean=exp_mean[tf_gene,1]
  RNA_mean_m=sum(RNA_mean)/length(RNA_mean)
  Motif_res1$Pat_mean_RNA[m]=RNA_mean_m

  RNA_ratio_m=length(which(aa[tf_gene,which(sub$AE=="Patient")]>0))/length(which(aa[tf_gene,which(sub$AE=="Patient")]>=0))##
  Motif_res1$Pat_RNA_ratio[m]=RNA_ratio_m
}
rownames(Motif_res1) <- Motif_res1$motif
write.csv(Motif_res1,paste0(cell[i],"_DEG_T-up-cisPeak_Motif_res.csv"),quote=F)

RES_all=rbind(RES_all,cbind(Motif_res1,cell[i]))
}
write.table(RES_all,"9cell_cis_Motif_enrichment_res.txt",sep="\t",quote=F)


##filter TF by enrich result and Motif expression
##Ranking Motifs by enrich pvalue
library(ggplot2)
library(ggrepel)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
#cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
Motif_res=read.table("9cell_cis_Motif_enrichment_res.txt",sep="\t",head=T)
##筛选chromvar_padj<0.01&chromvar_avg_diff>0.5
Motif_res1=Motif_res[Motif_res$enrich_padj<0.05,]

pdf("Fig4D_Motif_enrich_res1.pdf",height=6,width=5)
ggplot(Motif_res1,aes(y=motif.name,x=cell.i.,color=enrich_fold.enrichment,size=enrich_percent.observed))+
geom_point(alpha=0.8) +
labs(x= "Cell") +ylab("Motif") +
scale_color_gradient(low = "yellow", high = "darkblue")+
theme_bw()+
scale_size(breaks = c( 10,20, 30,40),range = c(0,6))
dev.off()


###再看一下这几个TF在对应细胞中的转录水平的情况
pdf("Fig4D_Motif-ranking_exp_FC.pdf",height=6,width=5)
p1 <- ggplot(Motif_res1,
  aes(y=motif.name,x=cell.i.,
  color=RNA_log2FC,size=Pat_RNA_ratio)) +
	geom_point(alpha=1) + #, size=2
	labs(x= "Cell") +ylab("Motif") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
  theme_bw()+
  scale_size(breaks = c(0,0.1,0.2,0.3,0.5,0.75),range = c(1,7))
print(p1)
dev.off()



5. Fig 4E GRN
#####找peak对应的Motif及BindSite
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
res_peak_up=cis[which(cis$threshold=="UP"),6]
length(res_peak_up)
##Peak list 输出成fa文件
peak1=gsub(":","-",res_peak_up)
## Plasma B markerPeak
Peak_bed=gsub("-","\t",peak1,perl=T)
write.table(Peak_bed,"9Tcell_DEG_cis_peaks.bed",sep="\t",row.names=F,col.names=F,quote=F)

#提取peak的DNA序列出来
#linux
#source activate Pyth3_R4
bedtools getfasta -fi /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-bed /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/9Tcell_DEG_cis_peaks.bed \
-fo /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/9Tcell_DEG_cis_peaks.fa

###计算Binding
#linux
mkdir /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/MOODS_res/

cd /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/633_Motifs_pfm/
for i in *_pwm.pfm
do
echo "$i"
moods-dna.py --m "$i" \
-s /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/9Tcell_DEG_cis_peaks.fa \
--p-value 0.0001  \
--outfile /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/MOODS_res/${i%.*}_motif_bind_0.0001.txt
done


###找之前的结果中比较关键的TF
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
res=read.table("9cell_cis_Motif_enrichment_res.txt",sep="\t",head=T)
res_mark1=res[which(res$enrich_fold.enrichment>1.5&res$enrich_padj<0.1
                    &res$chromvar_avg_diff>0&res$RNA_log2FC>0),]
key_TF=unique(res_mark1$motif.name)
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)


Motif=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/Motif_TF_name.txt",sep="\t",head=T)
key_TF=unique(res_mark1$motif.name)
key_TF_ID=Motif[Motif$motif.name%in%key_TF,2:3]
Bind_res3=c()
for(t in 1:length(key_TF_ID[,1])){
  TF_ID=key_TF_ID[t,1]
  TF_name=key_TF_ID[t,2]
  Bind_res=read.table(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/MOODS_res/",
                              TF_ID,"-",TF_name,"_pwm_motif_bind_0.0001.txt"),sep=",")
  Bind_res1=Bind_res[,c(1,3,4,5,6)]
  Bind_res2=cbind(Bind_res1,TF_ID,TF_name)
  colnames(Bind_res2)=c("PeakRanges","width","strand","BindingScore","BindingSeq","TF_ID","IF_name")
  Bind_res3=rbind(Bind_res3,Bind_res2)
}

dim(Bind_res3)
Bind_res4=Bind_res3[Bind_res3$BindingScore>0,]
dim(Bind_res4)

###找到每个Binding结果对应的peak和基因以及对应的celltype
Res1=c()
for(b in 1:length(Bind_res4[,1])){
  Res=cbind(Bind_res4[b,],cis[cis$Peak==Bind_res4[b,1],])
  Res1=rbind(Res1,Res)
}
Res1=Res1[Res1$BindingScore>7,]
Res2=unique(Res1[,c(7,8,9)])
write.table(Res1,"9cell_cis_Motif_enrichment_Binding_res_all.txt",quote=F,sep="\t",row.names=F)
write.table(Res2,"9cell_cis_Motif_enrichment_Binding_res_all_unique_last.txt",quote=F,sep="\t",row.names=F)

## network
#9cell_cis_Motif_enrichment_Binding_res_all.txt
res=read.table("9cell_cis_Motif_enrichment_Binding_res_all.txt",sep="\t",head=T)
res1=res[which(res$Celltype=="CD4 naive"&res$TF_name%in%c("CEBPG(var.2)","ATF4")),]
res2=res[which(res$Celltype=="CD4 Tm"&res$TF_name%in%c("MEF2A","NR3C1","NR3C2")),]
res3=res[which(res$Celltype=="CD8 Tcm"&res$TF_name%in%c("FOXK1","FOXK2","FOXN3","FOXO3","FOXP1","NR4A2","RELA")),]
res4=res[which(res$Celltype=="CD8 Tem"&res$TF_name%in%c("NR3C1","REL","RELA")),]

res_all=rbind(res1,res2,res3,res4)
dim(res_all)
write.table(res_all,"Fig4E_last_4cell_Motif_Bind_res.txt",sep="\t",quote=F)

###cytoscape构建network



6.Fig 4F network 中的gene的功能富集情况
##需要进一步整理
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/")
res=read.table("network_gene_function_res.txt",sep="\t",head=T)
res$Pvalue= -res$LogP
res$Ratio=res$InTerm_Count/res$Term_Count
res$Description=factor(res$Description,levels=res$Description)

library(ggplot2)
pdf("Fig4F_network_Function_enrichment_res.pdf",width=9,height=5)
ggplot(res,aes(x=Pvalue,y=Description,fill=InTerm_Count))+
geom_bar(stat="identity", width=.8, position = "dodge")+
labs(x="-Log10(P)",y="Term Description")+
theme_bw()
dev.off()














