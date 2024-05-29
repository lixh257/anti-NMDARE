FigS6

1.Fig S6A
library(Seurat)
library(ggplot2)
library(dittoSeq)
library(MySeuratWrappers)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
cell=table(combined_M$celltype2)
colnames(combined_M@meta.data)
table(combined_M@meta.data[,c(48,60)])
data=table(combined_M@meta.data[,c(48,60)])/apply(table(combined_M@meta.data[,c(48,60)]),1,sum)
data=as.data.frame(data)

pal_pat <- c("#756AB6","#AC87C5","#E0AED0")
pal_Nor = c("#A3D8F4","#8CC0DE","#A1CCD1","#64C9CF","#01A9B4","#176B87")
pdf("FigS6A_barplot_cellratio_of_Mcell_Sample.pdf",width=12,height=8)
ggplot(data,
aes(x=factor(celltype2,levels=c("CM","ncM","cDC","pDC")),
y=Freq,
fill=factor(Label,levels=c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6","Patient1","Patient3","Patient4"))))+
geom_bar(stat="identity", width=.8, position = "dodge")+
scale_fill_manual(values=c(pal_Nor,pal_pat),name="Sample")+
theme_bw()+
labs(x="Celltype",y="The proportion of each individual")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()


2. Fig S6B  DEG gene number
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig5/")
DEG_cM=read.csv("DEGs_CM_res.csv")
DEG_ncM=read.csv("DEGs_ncM_res.csv")
DEG_cDC=read.csv("DEGs_cDC_res.csv")
DEG_cM_UP=DEG_cM[which(DEG_cM$threshold=="UP"),1]
DEG_ncM_UP=DEG_ncM[which(DEG_ncM$threshold=="UP"),1]
DEG_cDC_UP=DEG_cDC[which(DEG_cDC$threshold=="UP"),1]
inter=intersect(DEG_cDC_UP,intersect(DEG_cM_UP,DEG_ncM_UP))

DEG_cM_Down=DEG_cM[which(DEG_cM$threshold=="Down"),1]
DEG_ncM_Down=DEG_ncM[which(DEG_ncM$threshold=="Down"),1]
DEG_cDC_Down=DEG_cDC[which(DEG_cDC$threshold=="Down"),1]
res=matrix(0,ncol=3,nrow=6)
res[1,]=c("CM","Up",length(DEG_cM_UP))
res[2,]=c("ncM","Up",length(DEG_ncM_UP))
res[3,]=c("cDC","Up",length(DEG_cDC_UP))
res[4,]=c("CM","Down",length(DEG_cM_Down))
res[5,]=c("ncM","Down",length(DEG_ncM_Down))
res[6,]=c("cDC","Down",length(DEG_cDC_Down))
colnames(res)=c("Celltype","Regu","number")
write.table(res,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/FigS6B_DEG_number_res.txt",sep="\t",quote=F,row.names=F)


res=read.table("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/FigS6B_DEG_number_res.txt",sep="\t",head=T)
library(ggplot2)
library(plyr)
library(gridExtra)
pal <- paletteer::paletteer_d("ggsci::category20_d3")
p_up=ggplot(subset(res,Regu=="Up"),aes(x=factor(Celltype,levels=c("CM","ncM","cDC")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))

p_down=ggplot(subset(res,Regu=="Down"),aes(x=factor(Celltype,levels=c(
"CM","ncM","cDC")),y=number,
fill=Celltype))+
geom_bar(position="stack",stat="identity")+
scale_fill_manual(values=pal,name="Cell Type")+
theme_bw()+
labs(x="Cell Type",y="Cell number")+
theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
## Plutting it together
pdf("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/FigS6B_DEG_number.pdf",width=12,height=5)
grid.arrange(p_up, p_down, ncol=2)
dev.off()



3-5.Fig S6C-E
#DEG function enrichment results
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/")
res_cM=read.table("cM_Function_enrichment_res.txt",sep="\t",head=T)
res_cM$Pvalue= -res_cM$LogP
res_cM$Ratio=res_cM$InTerm_Count/res_cM$Term_Count
res_cM$Description=factor(res_cM$Description,levels=res_cM$Description)

library(ggplot2)
pdf("FigS6C_cM_Function_enrichment_res.pdf",width=9,height=5)
ggplot(res_cM,aes(x=Pvalue,y=Description,fill=InTerm_Count))+
geom_bar(stat="identity", width=.8, position = "dodge")+
labs(x="-Log10(P)",y="Term Description")+
theme_bw()#+ scale_fill_fermenter( palette = "PuOr")#n.breaks = 9,
dev.off()


res_ncM=read.table("ncM_Function_enrichment_res.txt",sep="\t",head=T)
res_ncM$Pvalue= -res_ncM$LogP
res_ncM$Ratio=res_ncM$InTerm_Count/res_ncM$Term_Count
res_ncM$Description=factor(res_ncM$Description,levels=res_ncM$Description)

pdf("FigS6C_ncM_Function_enrichment_res.pdf",width=9,height=5)
ggplot(res_ncM,aes(x=Pvalue,y=Description,fill=InTerm_Count))+
geom_bar(stat="identity", width=.8, position = "dodge")+
labs(x="-Log10(P)",y="Term Description")+
theme_bw()#+ scale_fill_fermenter( palette = "PuOr")#n.breaks = 9,
dev.off()


res_cDC=read.table("cDC_Function_enrichment_res.txt",sep="\t",head=T)
res_cDC$Pvalue= -res_cDC$LogP
res_cDC$Ratio=res_cDC$InTerm_Count/res_cDC$Term_Count
res_cDC$Description=factor(res_cDC$Description,levels=res_cDC$Description)

pdf("FigS6C_cDC_Function_enrichment_res.pdf",width=9,height=5)
ggplot(res_cDC,aes(x=Pvalue,y=Description,fill=InTerm_Count))+
geom_bar(stat="identity", width=.8, position = "dodge")+
labs(x="-Log10(P)",y="Term Description")+
theme_bw()#+ scale_fill_fermenter( palette = "PuOr")#n.breaks = 9,
dev.off()




6. Fig S6F
####TLR 通路gene在Myeloid中的表达
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/FigS6/")
combined_M=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/2-Annotation/Merge_5lib_Myeloid.rds")
gene1=c("TIRAP","MYD88","IRAK4","IRAK1","IRAK2","TRAF6","TAB2","TAB3","MAP3K7","MAP2K3","MAP2K6",
        "MAP2K4","MAP2K7","MAPK14","MAPK8","NFKB1","NFKB2","REL","RELA",
        "CHUK")#"MKK3-MAP2K3","MKK6-MAP2K6","MKK4-MAP2K4","MKK7-MAP2K7",
        #p38-MAPK14,IKKa-CHUK,JNK-MAPK8,TAK1-MAP3K7

Idents(combined_M)="celltype2"
pdf("FigS6F_TLR_pathway_gene_RNA_Myeloid.pdf",width=5,height=7)
DotPlot(combined_M, features = gene1,dot.scale = 8)+
coord_flip()+#旋转图片  
theme_bw()+#去除背景，
theme(panel.grid = element_blank(),  
axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
#0099CCFF,#2171B5FF,#C6DBEFFF
scale_color_gradientn(values = seq(-1,1,0.2),colours =c("#253494FF","#2171B5FF","#C6DBEFFF",'#FFCC33'))+
#scale_color_gradientn(values = seq(-1,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
#颜色渐变设置  
scale_size(breaks = c(0,10, 20, 30, 50),range = c(-1,7))+
labs(y="Cell Type",x="IFNR Genes")+guides(size=guide_legend(order=3)) +
scale_y_discrete(limits =  c("CM","ncM","cDC")
)
dev.off()


DefaultAssay(combined_M)="RNA"
RNA_data=GetAssayData(combined_M,slot="data")
RNA_data1=RNA_data[gene1,]
cell=c("CM","ncM","cDC","pDC")
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
colnames(RES)=c("Gene","Cell","Pvalue","FC")
write.table(RES,"FigS6F_TLR_pathway_gene_Pat_vs_Nor_test_myeloid.txt",sep="\t",quote=F,row.names=F)

RES=read.table("FigS6F_TLR_pathway_gene_Pat_vs_Nor_test_myeloid.txt",sep="\t",head=T)
pdf("FigS6F_TLR_pathway_gene_test_myeloid.pdf",width=5,height=7)
ggplot(RES,aes(x=Cell,y=factor(Gene,levels=c("TIRAP","MYD88","IRAK4","IRAK1","IRAK2",
                          "TRAF6","TAB2","TAB3","MAP3K7","MAP2K3","MAP2K6",
                          "MAP2K4","MAP2K7","MAPK14","MAPK8","NFKB1","NFKB2",
                          "REL","RELA","CREB1","JUN","CHUK")),
  size=-log10(Pvalue),color=log2(FC)))+geom_point()+
	labs(x= "cell") +ylab("Gene") +
	scale_color_gradient2(low = "#1F4172", high = "#982176")+
  theme_bw()+
  scale_size(breaks = c(0,2,5,10,15),range = c(0,10))
dev.off()



7-8.Fig S6G-H CellChat
#构建cellchat project
library(CellChat)
library(patchwork)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig7/")
combined=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Merge_after_filter1.rds")
Normal=subset(x=combined,AE=="Normal")
AE=subset(x=combined,AE=="Patient")

data.input  <- Normal@assays$RNA@data
identity = data.frame(group =Normal$celltype2, row.names = names(Normal$celltype2)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat.Normal <- createCellChat(object = data.input)
cellchat.Normal <- addMeta(cellchat.Normal, meta = identity, meta.name = "labels")
cellchat.Normal <- setIdent(cellchat.Normal, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.Normal@idents)
groupSize <- as.numeric(table(cellchat.Normal@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
interaction_res=CellChatDB$interaction
write.table(interaction_res,"CellChat_interaction_res.txt",sep="\t",quote=F)


geneInfo_res=CellChatDB$geneInfo
write.table(geneInfo_res,"CellChat_geneInfo_res.txt",sep="\t",quote=F)

#CellChatDB$interaction[1:2,]
#head(CellChatDB$cofactor)
#head(CellChatDB$complex)
#head(CellChatDB$geneInfo)
#dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat.Normal@DB <- CellChatDB.use # set the used database in the object
#unique(CellChatDB$interaction$annotation)
#"Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact"
#################################预处理，首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。
cellchat.Normal <- subsetData(cellchat.Normal) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4) # do parallel
cellchat.Normal <- identifyOverExpressedGenes(cellchat.Normal)
cellchat.Normal <- identifyOverExpressedInteractions(cellchat.Normal)
cellchat.Normal <- projectData(cellchat.Normal, PPI.human)
#然后，我们通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信。
#Compute the communication probability and infer cellular communication network
cellchat.Normal<- computeCommunProb(cellchat.Normal)
cellchat.Normal <- filterCommunication(cellchat.Normal, min.cells = 10)
#Infer the cell-cell communication at a signaling pathway level
cellchat.Normal<- computeCommunProbPathway(cellchat.Normal)
#Calculate the aggregated cell-cell communication network
cellchat.Normal <- aggregateNet(cellchat.Normal)
####所有显示重要通信的信令路径都可以通过cellchat@netP$pathways访问。
#cellchat.P1@netP$pathways
saveRDS(cellchat.Normal, file = "cellchat_Normal.rds")


data.input  <- AE@assays$RNA@data
identity = data.frame(group =AE$celltype2, row.names = names(AE$celltype2)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat.AE<- createCellChat(object = data.input)
cellchat.AE <- addMeta(cellchat.AE, meta = identity, meta.name = "labels")
cellchat.AE <- setIdent(cellchat.AE, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.AE@idents)
groupSize <- as.numeric(table(cellchat.AE@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat.AE@DB <- CellChatDB.use # set the used database in the object
#unique(CellChatDB$interaction$annotation)
#"Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact"
#################################预处理，首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。
cellchat.AE<- subsetData(cellchat.AE) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4) # do parallel
cellchat.AE <- identifyOverExpressedGenes(cellchat.AE)
cellchat.AE <- identifyOverExpressedInteractions(cellchat.AE)
cellchat.AE <- projectData(cellchat.AE, PPI.human)
#然后，我们通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信。
#Compute the communication probability and infer cellular communication network
cellchat.AE<- computeCommunProb(cellchat.AE,raw.use=FALSE)
cellchat.AE <- filterCommunication(cellchat.AE, min.cells = 10)
#Infer the cell-cell communication at a signaling pathway level
cellchat.AE<- computeCommunProbPathway(cellchat.AE)
#Calculate the aggregated cell-cell communication network
cellchat.AE <- aggregateNet(cellchat.AE)
####所有显示重要通信的信令路径都可以通过cellchat@netP$pathways访问。
#cellchat.P1@netP$pathways
saveRDS(cellchat.AE, file = "cellchat_AE.rds")


setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig7/")
cellchat.Nor=readRDS("cellchat_Normal.rds")
cellchat.AE=readRDS("cellchat_AE.rds")
object.list <- list(Nor = cellchat.Nor, AE = cellchat.AE)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Part I: Predict general principles of cell-cell communication
##1.1
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("1-Merge_cellcaht_total_number.pdf",width=6, height=4)
gg1 + gg2
dev.off()

##1.2.
pdf("Merge_cellcaht_total_interaction_compare.pdf",width=12, height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
##1.3
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("Merge_cellcaht_total_interaction_compare_heatmap.pdf",width=12, height=6)
gg1 + gg2
dev.off()

##1.3.1.Differential number of interactions or interaction strength among different cell types
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("Merge_cellcaht_total_interaction_circor_number.pdf",width=15, height=7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 8, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("Merge_cellcaht_total_interaction_circor_strength.pdf",width=15, height=7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
   netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 15, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}
dev.off()


##1.4.Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
object.list[[i]]=netAnalysis_computeCentrality(
  object =object.list[[i]],
  slot.name = "netP",
  net = NULL,
  net.name = NULL
)
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("FigS6G_Merge_cellchat_total_interaction_strength_2D.pdf",width=10, height=4)
patchwork::wrap_plots(plots = gg)
dev.off()


#Part II: Identify the conserved and context-specific signaling pathways
#2.1
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("FigS6H_Merge_cellchat_pathway_signaling.pdf",width=8, height=6)
gg1 + gg2
dev.off()



#2.2
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10)
pdf("Merge_cellchat_pathway_signaling_outgoing.pdf",width=10, height=10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "GnBu")
pdf("Merge_cellchat_pathway_signaling_incoming.pdf",width=10, height=10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
pdf("Merge_cellchat_pathway_signaling_overall.pdf",width=10, height=10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()




