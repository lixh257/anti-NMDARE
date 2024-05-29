FigS5

library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(AnnotationDbi)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(stringr)
library(ggplot2)
library(ggrepel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS5")
combined_T=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_TNK.rds")
cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
cell=unique(cis$Celltype)

1.FigS5A DEG vocalnoplot
for(i in 1:length(cell)){
DEG=read.csv(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/DEGs_",cell[i],"_res.csv"))
#DEG$threshold1="Stable"
#DEG[which(DEG$threshold=="UP"&DEG$gene%in%cis[which(cis$Celltype==cell[i]),2]),9]="UP"
#DEG[which(DEG$threshold=="Down"&DEG$gene%in%cis[which(cis$Celltype==cell[i]),2]),9]="Down"
p=ggplot(DEG,
aes(x = log2FoldChange, y = -log10(padj), colour = threshold,label=gene)) +     
geom_point(alpha=0.8, size=0.8) +
scale_color_manual(values=c("#297CA0", "#DDDDDD","#EA5455"))+
#xlim(c(-1, 1)) +
ggtitle("Volcano plot of Patient vs Normal") +  
xlab("log2 fold change") +      
ylab("-log10 padj-value") +  
theme_bw()+   
#scale_y_continuous(limits = c(0,100)) +     
theme(legend.position = "right",           
plot.title = element_text(size = rel(1.5), hjust = 0.5),           
axis.title = element_text(size = rel(1.25))) 
DEG$label=ifelse(DEG$threshold!="Stable"&DEG$padj<1e-10,DEG$gene,"")
p+geom_text_repel(data = DEG, 
                  aes(x = log2FoldChange, 
                      y = -log10(padj), 
                      label = label),
                      max.overlaps=Inf,
                      size = 2,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"), 
                      segment.color = "#979797", 
                      show.legend = FALSE)                                                                
ggsave(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS5/FigS5A_",cell[i],"DEGs_Volcanoplot_log2FC_0.5_Term_2.pdf"),width=5,height=4)
}


2. Fig S5B upset of DEG gene
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS5/")
#cis=read.table("Tcell_9type_DEG_ciscor_DORC_Peak_res_200kb_cis_p0.01.txt",sep="\t",head=T)
cell=c("CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK") 
res5=c()
for(i in 1:length(cell)){  
DEG=read.csv(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/DEGs_",cell[i],"_res.csv"))
#order cell
res1=unique(DEG[which(DEG$threshold=="UP"),1])
res2=unique(DEG[which(DEG$threshold=="Down"),1])
res3=cbind(cell[i],"Up",res1)
res4=cbind(cell[i],"Down",res2)
res5=rbind(res5,res3,res4)
}
res6=res5[res5[,2]=="Up",]

cell=c("CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK") 
res=matrix(0,ncol=10,nrow=length(unique(res6[,3])))
res[,1]=unique(res6[,3])
colnames(res)=c("Gene",cell)
for(i in 1:length(res[,1])){
  cc=res6[which(res6[,3]==res[i,1]),1]
  res[i,cc]=1
}
write.table(res,"FigS5B_upset_res.txt",sep="\t",quote=F,row.names=F)

library(UpSetR)
res=read.table("FigS5B_upset_res.txt",sep="\t",head=T)
pdf("FigS5B_upset.pdf",width=10,height=8)
upset(res,nsets = 9)
dev.off()



3.Fig S5C
#####每个cell的motif ranking result
library(ggplot2)
library(ggrepel)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS5/")
motif_res=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig4/9cell_cis_Motif_enrichment_res.txt",sep="\t",head=T)
cell=c("CD4 naive","CD4 Tm","CD4 Treg","CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg","NK")  
for(i in 1:length(cell)){
    res1=motif_res[motif_res$cell.i.==cell[i],]
pdf(paste0(cell[i],"_Motif-ranking.pdf"),height=4,width=8)
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
}


