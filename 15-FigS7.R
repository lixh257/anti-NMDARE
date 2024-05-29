
FigS7
1.Fig S7A Model

2.Fig S7B
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS7/")
###TNF pathway gene expression
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(SeuratObject)
library(ggplot2)
library(devtools)
library(ggpubr)
library(dplyr)
Down_gene=c("TRADD","TRAF2",
       "MAP3K5","MAP3K4","MAP3K7","MAPK8","MAPK9","JUN",#AP-1
      "RIPK1","MAP3K3","MAP3K6", "MAPK1","MAPK3",#MAPK
      "MAP3K14","IKBKB","CHUK","IKBKB","RELA","REL","NFKB1","NFKB2"#NFKB
      )
      #IKKa-CHUK,IKKb-IKBKB,RIP-RIPK1,MAPK-MAPK1/MAPK3;
      #ASK1-MAP3K5;JNK-MAPK8/MAPK9/MAPK10,AP1-JUN
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
DefaultAssay(combined)="RNA"
RNA_data=GetAssayData(combined,slot="data")
RNA_data1=RNA_data[Down_gene,]

cell=c("CM","ncM","cDC","B naive","B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg",
        "CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg")
RES=c()
for(i in 1:length(Down_gene)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[Down_gene[i],which(combined$celltype2==cell[n]&combined$AE=="Patient")]
  exp_Nor=RNA_data1[Down_gene[i],which(combined$celltype2==cell[n]&combined$AE=="Normal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  exp_ratio=length(which(RNA_data1[Down_gene[i],which(combined$celltype2==cell[n])]>0))/length(RNA_data1[1,which(combined$celltype2==cell[n])])
  exp_mean=sum(RNA_data1[Down_gene[i],which(combined$celltype2==cell[n])])/length(RNA_data1[1,which(combined$celltype2==cell[n])])
  res=cbind(Down_gene[i],cell[n],Test_p,FC,exp_ratio,exp_mean)
  RES=rbind(RES,res)
}
}
colnames(RES)=c("Gene","Cell","Pvalue","FC","exp_ratio","exp_mean")
write.table(RES,"FigS7B_down_pathway_gene_Pat_vs_Nor_test_allcell.txt",sep="\t",quote=F,row.names=F)





RES=read.table("FigS7B_down_pathway_gene_Pat_vs_Nor_test_allcell.txt",sep="\t",head=T)
pdf("FigS7B_down_pathway_gene_Pat_vs_Nor_test_allcell.pdf",width=8,height=6)
ggplot(RES,aes(x=factor(Cell,levels=c("CM","ncM","cDC","B naive",
      "B memory","B plasma","CD4 naive","CD4 Tm","CD4 Treg",
       "CD8 naive","CD8 Tcm","CD8 Tem","CD8 MAIT","CD8 Treg")),y=factor(Gene,levels=c("TRADD","TRAF2",
       "MAP3K5","MAP3K4","MAP3K7","MAPK8","MAPK9","JUN",#AP-1
      "RIPK1","MAP3K3","MAP3K6", "MAPK1","MAPK3",#MAPK
      "MAP3K14","IKBKB","CHUK","RELA","REL","NFKB1","NFKB2"#NFKB
      )),
  size=exp_ratio,color=log2(FC)))+geom_point()+
	labs(x= "cell") +ylab("Gene") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
  theme_bw()+
  scale_size(breaks = c(0,0.1,0.2,0.3,0.4,0.5),range = c(0,10))
dev.off()