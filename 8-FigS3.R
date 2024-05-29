
Fig S3
################################################################################
1.Fig S3A 看一下GRN里面的几个TF的表达情况
#####
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS3/")
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
TF=c("IRF9","IRF4","IRF1","IRF7","IRF3","IRF2","STAT1::STAT2",
     "TCF12(var.2)","TWIST1","HSF4")#,"BHLHA15(var.2)"
Motif=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/Motif_TF_name.txt",sep="\t",head=T)
TF_gene=Motif[Motif$motif.name%in%TF,4]#12
DefaultAssay(BCell)="RNA"
Idents(BCell)="celltype2"
pdf("FigS3A_TF_gene_Bcell_exp.pdf",width=4,height=6)
DotPlot(BCell, features = TF_gene,dot.scale = 8)+
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


2.Fig S3B 几个TF具体的表达和chromvar活性的展示
##FigS3B_1
####Featureplot
library(dittoSeq)
TF_gene=c("IRF4","IRF1","IRF2","TCF12","STAT1","STAT2")
DefaultAssay(BCell)="SCT"###降低了差异幅度
sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
            "#2DB7A3","#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
for(i in 1:length(TF_gene)){
pdf(paste0("FigS3B_Bp_TF_Featureplot_",TF_gene[i],"_RNA.pdf"),width=4,height=3)
dittoDimPlot(BCell, 
             TF_gene[i],
             size = 1,
             min.color =sunrise[1:4],
             max.color =sunrise[5:9],
             show.axes.numbers = FALSE,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
}

##FigS3B_2
DefaultAssay(BCell) <- "SCT"
Idents(BCell)="celltype2"
cell=c("B naive","B memory","B plasma")
names(cell) <- levels(BCell)
BCell <- RenameIdents(BCell, cell)
pdf(paste0("FigS3B_TFgene_exp_mean_vlnplot_1.pdf"),width=6,height=4)
#par(mfrow = c(4, 4))
for(i in 1:length(TF_gene)){
p=VlnPlot(BCell, features = TF_gene[i], 
         pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols = paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20)]
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


##FigS3B_3
###chromvar
TF=c("IRF4","IRF1","IRF2","STAT1::STAT2","TCF12(var.2)")
TF_ID=unique(Motif[Motif$motif.name%in%TF,1])
DefaultAssay(BCell) <- "chromvar"
data1=GetAssayData(BCell)
Idents(BCell)="celltype2"
cell=c("B naive","B memory","B plasma")
names(cell) <- levels(BCell)
BCell <- RenameIdents(BCell, cell)
pdf(paste0("FigS3B_TF_Motif_chromvar_vlnplot.pdf"),width=4,height=3)
#par(mfrow = c(4, 4))
for(i in 1:length(TF_ID)){
p=VlnPlot(BCell, features = TF_ID[i], 
         pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols =paletteer::paletteer_d("ggsci::category20_d3")[c(10,1,20)]
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("CellType") + ylab("TF activity") + 
       #stat_compare_means( method = "wilcox.test") +
       #stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()

##FigS3B_4
solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D')
DefaultAssay(BCell)="chromvar"
for(i in 1:length(TF_ID)){
pdf(paste0("FigS3B_Bp_TF_Featureplot_",TF_ID[i],"_chromvar.pdf"),width=4,height=3)
dittoDimPlot(BCell, 
             TF_ID[i],
             size = 1,
             min.color =solarExtra[1:5],
             max.color = solarExtra[6:9],
             show.axes.numbers = FALSE,
             order = "increasing",#order = c("unordered", "increasing", "decreasing", "randomize")
             reduction = "wnn.umap",
             labels.size = 50
             )                
dev.off()
}

##FigS3B_5
##展示Motif logo
DefaultAssay(BCell)="chromvar"
#a=GetAssayData(BCell)
for(i in 1:5){
pdf(paste0("FigS3B_plotmotif_",TF_ID[i],"_5.pdf"),width=5,height=2)
MotifPlot(
  object = BCell,
  motifs = TF_ID[i],
  assay = 'Peaks1'
)
dev.off()
}



3. Fig S3C Some gene expression in Plasma B cells

gene=c("IRF1","IRF2","IRF4","TCF12","STAT1","STAT2","IFNAR1","IFNAR2","TYK2","JAK1")
Bp=subset(BCell,celltype2=="B plasma")
DefaultAssay(Bp)="RNA"
Idents(Bp)="Label"
levels(Bp)=c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Patient1","Patient3","Patient4")

pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#539092","#8BCFCC","#51ADCF","#0278AE","#01A9B4","#A3D8F4")

pdf(paste0("FigS3C_Bplasma_TF_gene_PvsN.pdf"),width=6,height=3)
for(i in 1:length(gene)){
p=VlnPlot(Bp, features = gene[i], 
         pt.size = 0.1, 
        combine = T
        ) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("Samples") + ylab("Gene Expression") + 
       scale_fill_manual(values=c(pal_Nor,pal_pat),name="Samples")+
       theme(text = element_text(size = 10))+ #, angle = 45
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
print(p) 
}
dev.off()

#计算P值：
DefaultAssay(Bp)="RNA"
RNA_data=GetAssayData(Bp,slot="data")
result1=c()
for(i in 1:length(gene)){
  TT=wilcox.test(RNA_data[gene[i],which(Bp$AE=="Patient")],RNA_data[gene[i],which(Bp$AE=="Normal")])
  TTP=TT[[3]]
  result=cbind(gene[i],TTP)
  result1=rbind(result1,result)
}


4.Fig S3D
#######################几个TF是由什么TF/位点调控的呢？
TF=c("IRF1","IRF4","IRF2","STAT1","STAT2","TCF12")
cis=read.table("Bcell_all_cisCor_filter_200kb.txt",sep="\t",head=T)
peak=cis[which(cis$Gene%in%TF),2]
##Peak list 输出成fa文件
peak1=gsub(":","-",peak)
## Plasma B markerPeak
Peak_bed=gsub("-","\t",peak1,perl=T)
write.table(Peak_bed,"IRF4_ex_TF_cis_peaks.bed",sep="\t",row.names=F,col.names=F,quote=F)

#提取peak的DNA序列出来
#linux
source activate Pyth3_R4
bedtools getfasta -fi /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-bed /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_TF_cis_peaks.bed \
-fo /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_TF_cis_peaks.fa

###计算Binding
#linux
mkdir /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_IRF4_ex_TF_cis_peaks/
cd /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/633_Motifs_pfm/
for i in *_pwm.pfm
do
echo "$i"
moods-dna.py --m "$i" \
-s /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_TF_cis_peaks.fa \
--p-value 0.001  \
--outfile /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_IRF4_ex_TF_cis_peaks/${i%.*}_motif_bind_0.001.txt
done


###
DEG=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/DEG_celltype2_Bcell_padj0.01_log2FC0.5.txt",sep="\t",head=T)
Moitf=read.table("Motif_TF_name.txt",sep="\t",head=T)
TF=Moitf[,4]
TF_gene=intersect(TF,DEG[DEG$cluster=="B plasma",7])
# [1] "JUN"   "RUNX2" "MEF2C" "ESR1"  "TP63"  "CUX1"  "ELF1"  "PRDM1" "TCF12"
#[10] "IRF1"  "FOXN3" "XBP1"  "IRF4"
###如果不要求TF必须是差异表达的
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/")
Moitf=read.table("Motif_TF_name.txt",sep="\t",head=T)
Motif1=unique(Motif[,c(1,3)])
Bind_res3=c()
for(n in c(1:11,13:63,65:81,83:98,100:126,128:137,139,141,143:155,158:270,
           272:279,281:285,287:345,347:359,361:424,426:429,431:432,
           434:466,468:488,490:496,498:513,515,517:593,596:length(Motif1[,1]))){##有些没有找到结合位点所以文件是0
  print(n)
  TF_ID=Motif1[n,1]
  TF_name=Motif1[n,2]
  Bind_res=read.table(paste0("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/MOODS_res_IRF4_ex_TF_cis_peaks/",
                              TF_ID,"-",TF_name,"_pwm_motif_bind_0.001.txt"),sep=",")
  Bind_res1=Bind_res[,c(1,3,4,5,6)]
  Bind_res2=cbind(Bind_res1,TF_ID,TF_name)
  colnames(Bind_res2)=c("PeakRanges","width","strand","BindingScore","BindingSeq","TF_ID","IF_name")
  Bind_res3=rbind(Bind_res3,Bind_res2)
}
write.table(Bind_res3,"IRF4_ex_Binding_site_res.txt",sep="\t",quote=F,row.names=F)
#cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_Binding_site_res.txt /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS3
cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_TF_cis_peaks.fa /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS3
cp -r /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig2/IRF4_ex_TF_cis_peaks.bed /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS3
###通过Binding score以及这些TF在BCell中的表达情况，最终发现IRF4的调控TF可能是TCF12，IRF1和EWSR1-FLI1

###看一下track的情况
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS3/")
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
gene=c("FNDC3B","IRF1","IRF2","STAT1","STAT2","IRF4","TCF12")
BCell=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Bcell.rds")
DefaultAssay(BCell)="Peaks1"
BCell <- RegionStats(BCell, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(BCell) <- 'SCT' 
SCT=GetAssayData(BCell,slot="data")
Gene=gene[which(gene%in%rownames(SCT))]
BCell <- LinkPeaks(
  object = BCell,
  peak.assay = "Peaks1",
  expression.assay = "SCT",
  genes.use = Gene,
  pvalue_cutoff = 0.01,
  score_cutoff = 0
  )
Idents(BCell) <- "celltype2"
levels(BCell) <- c("B naive","B memory","B plasma")
DefaultAssay(BCell)="Peaks1"

Bind=read.table("IRF4_ex_Binding_site_res.txt",sep="\t",head=T)
cis=read.table("Bcell_all_cisCor_filter_200kb.txt",sep="\t",head=T)
peak1=cis[which(cis$Gene=="IRF4"),2]

Sites=peak1
Sites=gsub(":","-",Sites)
ranges.show = StringToGRanges(Sites)
pdf(paste0("FigS3D_TF_gene_track_200kb_IRF4_new.pdf"),height=6,width=10)
p=CoveragePlot(
   object = BCell,
   region = "IRF4", 
   region.highlight =ranges.show,
   #ranges = Peaks1,
   features = "IRF4",
   expression.assay = "SCT",
   extend.upstream = 180000,
   extend.downstream = 1000,
   #ymax = 20,##change normalized signal range
   links=T,
   window = 200
 )
 print(p)
 dev.off()



###展示seqence:在fa文件里面
##展示Motif logo
DefaultAssay(BCell)="chromvar"
#TCF12-MA1648.1
#EWSR1-FLI1-MA0149.1
#IRF1-MA0050.2
#a=GetAssayData(BCell)
TF_ID=c("MA1648.1","MA0149.1","MA0050.2")
TF=c("TCF12","EWSR1-FLI1","IRF1")
for(i in 1:3){
pdf(paste0("FigS3D_plotmotif_",TF[i],"_new.pdf"),width=5,height=3)
MotifPlot(
  object = BCell,
  motifs = TF_ID[i],
  assay = 'Peaks1'
)
dev.off()
}



