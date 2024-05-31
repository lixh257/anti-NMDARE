library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(edgeR)
library(Matrix)
library(ggplot2)
library(scater)
library(scuttle)
library(DropletUtils)
library(scDblFinder)
library(BiocParallel)
library(RColorBrewer)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(eulerr)
library(pheatmap)
library(parallel)
library(tidyverse)
library(vcfR)
library(pheatmap)
library(ggbeeswarm)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", "#9e4747", 
	"#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8", "#9AD3DD", "#3071A9", "#31B3C2")
text_size <- 11
title_size <- 13
tag_thm <- theme(plot.tag=element_text(size=20, colour="black"), plot.margin=margin())

sce_patient <- readRDS("patient.rds")
sample_cls <- read.delim("clusters.tsv")
sample_cls <- sample_cls[which(sample_cls$status == "singlet"),]
sample_cls$assignment[which(sample_cls$assignment == "0")] <- "Patient5"
sample_cls$assignment[which(sample_cls$assignment == "1")] <- "Patient4"
sample_cls$assignment[which(sample_cls$assignment == "2")] <- "Patient1"
sample_cls$assignment[which(sample_cls$assignment == "3")] <- "Patient2"
sample_cls$assignment[which(sample_cls$assignment == "4")] <- "Patient3"
sample_cls <- sample_cls[which(sample_cls$assignment != "Patient5"),]
sce_patient <- subset(sce_patient, cells=intersect(sample_cls$barcode, colnames(sce_patient)))
celltype <- c("BNaive", "BMemory", "BPlasma", "Mono", "DC", "MKC", "TC")
cell_marker <- read.delim("cell_marker.txt")
sce_patient$cell.type <- factor(sce_patient$cell.type, levels=celltype)

sce_normal <- readRDS("normal.rds")
celltype <- c("BNaive", "BMemory", "BPlasma", "Mono", "DC", "MKC", "TC")
cell_marker <- read.delim("cell_marker.txt")
sce_normal$cell.type <- factor(sce_normal$cell.type, levels=celltype)

paa <- DotPlot(sce_patient, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"))
pab <- UMAPPlot(sce_patient, group.by="cell.type", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour=NULL)

pba <- DotPlot(sce_normal, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"))
pbb <- UMAPPlot(sce_normal, group.by="cell.type", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour=NULL)


colors <- colorRampPalette(brewer.pal(n=9, name="Reds"))(100)
split_mat <- readRDS("cls_rna_mat.rds")
ids <- which(split_mat[["info"]]$Group != "0")
split_mat[["mm"]] <- split_mat[["mm"]][, ids]
split_mat[["anno"]][["Group"]] <- split_mat[["anno"]][["Group"]][-1]
split_mat[["info"]] <- split_mat[["info"]][ids,, drop=F]
pca <- pheatmap(split_mat[["mm"]], cluster_cols=F, cluster_rows=F, show_rownames=F, show_colnames=F, 
	color=colors, annotation_colors=split_mat[["anno"]], border_color=NA, 
	labels_row="Mutation", labels_col="Cell", annotation_col=split_mat[["info"]])
split_cor <- readRDS("cls_cor_mat.rds")
ids <- which(split_cor[["info"]]$Group != "0")
split_cor[["mm"]] <- split_cor[["mm"]][ids, ids]
split_cor[["anno"]][["Group"]] <- split_cor[["anno"]][["Group"]][-1]
split_cor[["info"]] <- split_cor[["info"]][ids,, drop=F]
pcb <- pheatmap(split_cor[["mm"]], cluster_cols=F, cluster_rows=F, show_rownames=F, show_colnames=F, 
	color=colors, annotation_colors=split_cor[["anno"]], border_color=NA, 
	labels_row="Cell", labels_col="Cell", annotation_col=split_cor[["info"]], annotation_row=split_cor[["info"]])
pc <- grid.arrange(grobs=list(pca[[4]], pcb[[4]]), nrow=1, width=c(2, 1))

pheatmap(split_mat[["mm"]], cluster_cols=F, cluster_rows=F, show_rownames=F, show_colnames=F, 
	color=colors, annotation_colors=split_mat[["anno"]], cellheight=2, cellwidth=0.5, fontsize=10, border_color=NA, 
	labels_row="Mutation", labels_col="Cell", annotation_col=split_mat[["info"]], filename = "cls_rna.pdf")

pheatmap(split_cor[["mm"]], cluster_cols=F, cluster_rows=F, show_rownames=F, show_colnames=F, 
	color=colors, annotation_colors=split_cor[["anno"]], cellheight=0.3, cellwidth=0.5, fontsize=10, border_color=NA, 
	labels_row="Cell", labels_col="Cell", annotation_col=split_cor[["info"]], 
	annotation_row=split_cor[["info"]], filename = "cls_rna_cor.pdf")


mm_atac <- readRDS("mm_atac.rds")
mm_rna <- readRDS("mm_rna.rds")
mm_test <- readRDS("mm_test.rds")
mm <- data.frame(scale(mm_atac))
ca <- data.frame()
for (i in 1:nrow(mm)) ca <- rbind(ca, data.frame(X=rownames(mm)[i], Y=colnames(mm), Val=as.numeric(mm[i,])))
#ca$X <- factor(as.character(ca$X), levels=as.character(c(2,3,4,1,0)))
ca <- ca[which(ca$X != "0" & ca$Y != "Paient5"),]
pda <- ggplot(ca, aes(x=X, y=Y, fill=Val))+geom_tile()+theme_minimal()+
	labs(title=NULL, x="5'RNA", y="ATAC", fill="Inter\n(scaled)")+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	scale_fill_gradient2(low=brewer.pal(11,"Spectral")[11], mid=brewer.pal(11,"Spectral")[6], 
	high=brewer.pal(11,"Spectral")[1])+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks=element_blank(), panel.background=element_blank(), legend.key.height=unit(40, "pt"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"))

mm <- data.frame(scale(mm_rna))
cb <- data.frame()
for (i in 1:nrow(mm)) cb <- rbind(cb, data.frame(X=rownames(mm)[i], Y=colnames(mm), Val=as.numeric(mm[i,])))
#cb$X <- factor(as.character(cb$X), levels=as.character(c(2,3,4,1,0)))
cb <- cb[which(cb$X != "0" & cb$Y != "Paient5"),]
pdb <- ggplot(cb, aes(x=X, y=Y, fill=Val))+geom_tile()+theme_minimal()+
	labs(title=NULL, x="5'RNA", y="3'RNA", fill="Inter\n(scaled)")+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	scale_fill_gradient2(low=brewer.pal(11,"Spectral")[11], mid=brewer.pal(11,"Spectral")[6], 
	high=brewer.pal(11,"Spectral")[1])+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks=element_blank(), panel.background=element_blank(), legend.key.height=unit(40, "pt"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"))

mm <- data.frame(scale(t(mm_test)))
cc <- data.frame()
for (i in 1:nrow(mm)) cc <- rbind(cc, data.frame(X=rownames(mm)[i], Y=colnames(mm), Val=as.numeric(mm[i,])))
cc <- cc[which(cc$X != "Paient5" & cc$Y != "Paient5"),]
pdc <- ggplot(cc, aes(x=X, y=Y, fill=Val))+geom_tile()+theme_minimal()+
	labs(title=NULL, x="3'RNA", y="ATAC", fill="Inter\n(scaled)")+
	scale_x_discrete(expand=c(0, 0))+scale_y_discrete(expand=c(0, 0))+
	scale_fill_gradient2(low=brewer.pal(11,"Spectral")[11], mid=brewer.pal(11,"Spectral")[6], 
	high=brewer.pal(11,"Spectral")[1])+
	theme(axis.text=element_text(size=text_size, colour="black"), axis.title=element_text(size=title_size, colour="black"), 
	axis.ticks=element_blank(), panel.background=element_blank(), legend.key.height=unit(40, "pt"), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black"))

sce_bcell <-subset(sce_patient, cells=colnames(sce_patient)[grep("^B", sce_patient$cell.type)])
sce_bcell <- subset(sce_bcell, cells=intersect(sample_cls$barcode, colnames(sce_bcell)))
sce_bcell$cell.sample <- sample_cls$assignment[match(colnames(sce_bcell), sample_cls$barcode)]
pe <- wrap_elements(UMAPPlot(sce_bcell, group.by="cell.sample", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour=NULL))+tag_thm

sce_normal <- readRDS("normal.rds")
bcrinfo <- read.csv("data/normal_bcr.csv", h=T)
celltype <- c("BNaive", "BMemory", "BPlasma")
bcrinfo$type <- ""
for (type in celltype)
{
	bcs <- colnames(sce_normal)[which(sce_normal$cell.type == type)]
	for (bc in bcs)
	{
		ids <- which(bcrinfo$barcode == bc)
		if (length(ids) > 0) bcrinfo$type[ids] <- type
	}
}
bcrinfo <- bcrinfo[which(bcrinfo$type != ""),]
bcrinfo$raw_clonotype_id <- paste0(bcrinfo$type, "+", bcrinfo$raw_clonotype_id)
bcrinfo$raw_clonotype_id <- gsub("\\..*,", "_", paste0(bcrinfo$barcode, ",", bcrinfo$raw_clonotype_id))
bcrinfo$cc <- 0
clcs <- data.frame(table(gsub("^.*,", "", unique(paste0(bcrinfo$barcode, ",", bcrinfo$raw_clonotype_id)))))
for (i in 1:nrow(clcs)) bcrinfo$cc[which(bcrinfo$raw_clonotype_id == clcs[i, 1])] <- clcs[i, 2]
#bcrinfo <- bcrinfo[which(bcrinfo$c_gene != "" & bcrinfo$c_gene != "IGKC" & bcrinfo$c_gene != "None"),]
#if (length(grep("IGLC", bcrinfo$c_gene)) > 0) bcrinfo <- bcrinfo[-grep("IGLC", bcrinfo$c_gene),]
bcrinfo_normal <- bcrinfo
sce_patient <- readRDS("patient.rds")
sample_cls <- read.delim("clusters.tsv")
sample_cls <- sample_cls[which(sample_cls$status == "singlet"),]
sample_cls$assignment[which(sample_cls$assignment == "0")] <- "Patient5"
sample_cls$assignment[which(sample_cls$assignment == "1")] <- "Patient4"
sample_cls$assignment[which(sample_cls$assignment == "2")] <- "Patient1"
sample_cls$assignment[which(sample_cls$assignment == "3")] <- "Patient2"
sample_cls$assignment[which(sample_cls$assignment == "4")] <- "Patient3"
sample_cls <- sample_cls[which(sample_cls$assignment != "Patient5"),]
sce_patient <- subset(sce_patient, cells=intersect(sample_cls$barcode, colnames(sce_patient)))
sce_patient$cell.sample <- sample_cls$assignment[match(colnames(sce_patient), sample_cls$barcode)]
bcrinfo <- read.csv("data/patient_bcr.csv", h=T)
celltype <- c("BNaive", "BMemory", "BPlasma")
bcrinfo$type <- ""
for (type in celltype)
{
	bcs <- colnames(sce_patient)[which(sce_patient$cell.type == type)]
	for (bc in bcs)
	{
		ids <- which(bcrinfo$barcode == bc)
		if (length(ids) > 0) bcrinfo$type[ids] <- type
	}
}
cellsample <- unique(sce_patient$cell.sample)
bcrinfo$sample <- ""
for (sample in cellsample)
{
	bcs <- colnames(sce_patient)[which(sce_patient$cell.sample == sample)]
	for (bc in bcs)
	{
		ids <- which(bcrinfo$barcode == bc)
		if (length(ids) > 0) bcrinfo$sample[ids] <- sample
	}
}
bcrinfo <- bcrinfo[which(bcrinfo$type != "" | bcrinfo$sample != ""),]
bcrinfo$raw_clonotype_id <- paste0(bcrinfo$sample, "+", bcrinfo$type, "+", bcrinfo$raw_clonotype_id)
bcrinfo$raw_clonotype_id <- gsub("^.*,", "", paste0(bcrinfo$barcode, ",", bcrinfo$raw_clonotype_id))
bcrinfo$cc <- 0
clcs <- data.frame(table(gsub("^.*,", "", unique(paste0(bcrinfo$barcode, ",", bcrinfo$raw_clonotype_id)))))
for (i in 1:nrow(clcs)) bcrinfo$cc[which(bcrinfo$raw_clonotype_id == clcs[i, 1])] <- clcs[i, 2]
#bcrinfo <- bcrinfo[which(bcrinfo$c_gene != "" & bcrinfo$c_gene != "IGKC" & bcrinfo$c_gene != "None"),]
#if (length(grep("IGLC", bcrinfo$c_gene)) > 0) bcrinfo <- bcrinfo[-grep("IGLC", bcrinfo$c_gene),]
bcrinfo_patient <- bcrinfo

sce <- sce_normal
bcrinfo <- bcrinfo_normal
celltype <- c("BMemory")
res_total <- data.frame()
bcs <- c()
terms <- c()
for (type in celltype) bcs <- c(bcs, colnames(sce)[which(sce$cell.type == type)])
for (bc in bcs)
{
	ids <- which(bcrinfo$barcode == bc)
	if (length(ids) > 0) terms <- c(terms, ids)
}
bcrinfo <- bcrinfo[terms,]
res <- unique(bcrinfo$raw_clonotype_id)
res <- data.frame(Var1=res, Freq=bcrinfo$cc[match(res, bcrinfo$raw_clonotype_id)])
res$group <- "Unique"
res$group[which(res$Freq >= 2)] <- "2-4"
res$group[which(res$Freq >= 5)] <- "5-9"
res$group[which(res$Freq >= 10)] <- "≥10"
res$group <- factor(res$group, levels=c("Unique", "2-4", "5-9", "≥10"))
res <- data.frame(table(res$group))
res$Var1 <- factor(res$Var1, levels=c("Unique", "2-4", "5-9", "≥10"))
pea <- ggplot(res, aes(x=3, y=Freq, fill=Var1))+geom_col()+labs(title="BMemory\n(Normal)", x=NULL, y=NULL, fill="Num. of clones")+
	geom_text(aes(label=Freq), position=position_stack(vjust=0.5))+
	coord_polar(theta="y")+xlim(c(0.2, 3.5))+scale_fill_manual(values=col_list)+
	theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
	axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
	axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
	legend.text=element_text(colour="black", size=12))
sce <- sce_patient
sce$cell.temp <- paste0(sce$cell.sample, "+", sce$cell.type)
bcrinfo <- bcrinfo_patient
samples <- unique(bcrinfo$sample)
samples <- samples[order(samples)]
plsa <- lapply(samples, function(sample) 
{
	type <- paste0(sample, "+BMemory")
	res_total <- data.frame()
	bcs <- c()
	terms <- c()
	bcs <- c(bcs, colnames(sce)[which(sce$cell.temp == type)])
	for (bc in bcs)
	{
		ids <- which(bcrinfo$barcode == bc)
		if (length(ids) > 0) terms <- c(terms, ids)
	}
	bcrinfo <- bcrinfo[terms,]
	res <- unique(bcrinfo$raw_clonotype_id)
	res <- data.frame(Var1=res, Freq=bcrinfo$cc[match(res, bcrinfo$raw_clonotype_id)])
	res$group <- "Unique"
	res$group[which(res$Freq >= 2)] <- "2-4"
	res$group[which(res$Freq >= 5)] <- "5-9"
	res$group[which(res$Freq >= 10)] <- "≥10"
	res$group <- factor(res$group, levels=c("Unique", "2-4", "5-9", "≥10"))
	res <- data.frame(table(res$group))
	res$Var1 <- factor(res$Var1, levels=c("Unique", "2-4", "5-9", "≥10"))
	return(ggplot(res, aes(x=3, y=Freq, fill=Var1))+geom_col()+
		labs(title=paste0("BMemory\n(", sample, ")"), x=NULL, y=NULL, fill="Num. of clones")+
		geom_text(aes(label=Freq), position=position_stack(vjust=0.5))+
		coord_polar(theta="y")+xlim(c(0.2, 3.5))+scale_fill_manual(values=col_list)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12)))
})
sce <- sce_patient
sce$cell.temp <- paste0(sce$cell.sample, "+", sce$cell.type)
bcrinfo <- bcrinfo_patient
samples <- unique(bcrinfo$sample)
samples <- samples[order(samples)]
plsb <- lapply(samples, function(sample) 
{
	type <- paste0(sample, "+BPlasma")
	res_total <- data.frame()
	bcs <- c()
	terms <- c()
	bcs <- c(bcs, colnames(sce)[which(sce$cell.temp == type)])
	for (bc in bcs)
	{
		ids <- which(bcrinfo$barcode == bc)
		if (length(ids) > 0) terms <- c(terms, ids)
	}
	bcrinfo <- bcrinfo[terms,]
	res <- unique(bcrinfo$raw_clonotype_id)
	res <- data.frame(Var1=res, Freq=bcrinfo$cc[match(res, bcrinfo$raw_clonotype_id)])
	res$group <- "Unique"
	res$group[which(res$Freq >= 2)] <- "2-4"
	res$group[which(res$Freq >= 5)] <- "5-9"
	res$group[which(res$Freq >= 10)] <- "≥10"
	res$group <- factor(res$group, levels=c("Unique", "2-4", "5-9", "≥10"))
	res <- data.frame(table(res$group))
	res$Var1 <- factor(res$Var1, levels=c("Unique", "2-4", "5-9", "≥10"))
	return(ggplot(res, aes(x=3, y=Freq, fill=Var1))+geom_col()+
		labs(title=paste0("BPlasma\n(", sample, ")"), x=NULL, y=NULL, fill="Num. of clones")+
		geom_text(aes(label=Freq), position=position_stack(vjust=0.5))+
		coord_polar(theta="y")+xlim(c(0.2, 3.5))+scale_fill_manual(values=col_list)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12)))
})
pf <- wrap_elements(wrap_plots(c(plsa, list(pea), plsb), ncol=5)+plot_layout(guides="collect"))+tag_thm

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(wrap_elements(paa)+tag_thm, wrap_elements(pab)+tag_thm, wrap_elements(pba)+tag_thm, 
		wrap_elements(pbb)+tag_thm), nrow=1, widths=c(2,3.5,1.5,3.5))+plot_annotation(tag_levels=list(c("A", "", "B", "")))), 
	wrap_elements(wrap_elements(wrap_plots(pc))+plot_annotation(tag_levels=list(c("C")))), 
	wrap_elements(wrap_plots(list(pda, pdb, pdc), nrow=1)+tag_thm+plot_annotation(tag_levels=list(c("D", "", "", "")))), 
	wrap_elements(wrap_plots(list(pe, pf), nrow=1, widths=c(1,2.5))+tag_thm+plot_annotation(tag_levels=list(c("E", "F"))))
	), ncol=1, heights=c(2,2,2,2)), 
	width=20, height=21, dpi=200, filename="bcr_figs.png")

ggsave(plot=wrap_plots(list(
	wrap_elements(wrap_plots(list(wrap_elements(paa)+tag_thm, wrap_elements(pab)+tag_thm, wrap_elements(pba)+tag_thm, 
		wrap_elements(pbb)+tag_thm), nrow=1, widths=c(2,3.5,1.5,3.5))+plot_annotation(tag_levels=list(c("A", "", "B", "")))), 
	wrap_elements(wrap_elements(wrap_plots(pc))+plot_annotation(tag_levels=list(c("C")))), 
	wrap_elements(wrap_plots(list(pda, pdb, pdc), nrow=1)+tag_thm+plot_annotation(tag_levels=list(c("D", "", "", "")))), 
	wrap_elements(wrap_plots(list(pe, pf), nrow=1, widths=c(1,2.5))+tag_thm+plot_annotation(tag_levels=list(c("E", "F"))))
	), ncol=1, heights=c(2,2,2,2)), 
	width=20, height=21, dpi=200, filename="bcr_figs.pdf")

