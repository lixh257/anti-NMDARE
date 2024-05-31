##BCR-seq--Patient
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
library(pheatmap)
library(ggheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(ggbeeswarm)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
text_size <- 11
title_size <- 13
tag_thm <- theme(plot.tag=element_text(size=14, colour="black", face="bold"), plot.margin=margin())


#fileList <- list.files(path="data", pattern="*.h5$")
#sceList <- lapply(fileList, function(x){CreateSeuratObject(Read10X_h5(paste0("data/", x)), project=gsub("_rna.h5", "", x))})
#sce <- merge(sceList[1][[1]], y=sceList[-1], add.cell.ids=gsub(".h5", "", fileList), project="striatum_RNA")
sce <- CreateSeuratObject(Read10X_h5("data/patient_rna.h5"), project="patient")
sce[["cell.sample"]] <- Idents(sce)
sce[["cell.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
cs_info <- data.frame(sce@meta.data)
pa <- ggplot(cs_info, aes(x="Term", y=nFeature_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Genes", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pb <- ggplot(cs_info, aes(x="Term", y=nCount_RNA))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Counts", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pc <- ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
ggsave(plot=pa+pb+pc, width=6, height=4, dpi=200, "rna_qc_norm.png")

pa <- ggplot(cs_info, aes(x="Term", y=log10(nFeature_RNA)))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+scale_y_continuous(breaks=seq(2, 4, 0.5), limits=c(2, 4))+
	labs(title=NULL, x="Genes (log10)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pb <- ggplot(cs_info, aes(x="Term", y=log10(nCount_RNA)))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Counts (log10)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pc <- ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_y_continuous(breaks=seq(5, 25, 5), limit=c(0, 30))+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
ggsave(plot=pa+pb+pc, width=6, height=4, dpi=200, "rna_qc.png")
ggsave(plot=pa+pb+pc, width=6, height=4, dpi=200, "rna_qc.pdf")


sce <- as.SingleCellExperiment(sce)
id <- grep("^IG[LHK].*", rownames(sce))
if (length(id) > 0) sce <- sce[-id,]
sce <- scDblFinder(sce, samples="cell.sample")
#sce$scDblFinder.class <- "singlet"
dbi_info <- sce$scDblFinder.class
mt_limit <- 10
mt_info <- data.frame(value=sce$cell.mt, type="F")
mt_info$type[which(mt_info$value < mt_limit)] <- "P"
cs_info <- data.frame(value=colSums(counts(sce) > 0), umi=colSums(counts(sce)), type="F")
cs_limit <- c(700, 3500)
cs_info$type[which(cs_info$value > cs_limit[1] & cs_info$umi > cs_limit[1]*2 & cs_info$value < cs_limit[2] & cs_info$umi < cs_limit[2]*3)] <- "P"
sce <- sce[, which(mt_info$type == "P" & cs_info$type == "P" & sce$scDblFinder.class == "singlet")]
sce <- sce[rowSums(counts(sce) > 0) > 10,]
pa <- ggplot(mt_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_y_continuous(breaks=seq(10, 90, 10))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=mt_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pb <- ggplot(cs_info, aes(x="Term", y=umi, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Counts", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1]*2, seq(5000, 50000, 5000)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=c(cs_limit[1]*2, cs_limit[2]*3))+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pc <- ggplot(cs_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Feature", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1], seq(500, 5000, 500)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=c(cs_limit[1], cs_limit[2]))+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())

sce <- addPerCellQC(sce, percent_top=c(20, 50, 100, 200))
sce$total_features <- sce$detected
sce$log10_total_features <- log10(sce$detected)
sce$total_counts <- sce$sum
sce$log10_total_counts <- log10(sce$sum)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features
mod <- loess(colData(sce)$log10_total_features~colData(sce)$log10_total_counts)
pred <- predict(mod, newdata=data.frame(log10_total_counts=colData(sce)$log10_total_counts))
sce$featcount_dist <- colData(sce)$log10_total_features - pred
sce$pct_counts_top_20_features <- colData(sce)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(sce)))[[1]]]]
sce$pct_counts_top_50_features <- colData(sce)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(sce)))[[1]]]]
vars <- c("log10_total_counts:both:5", "log10_total_features:both:5", "pct_counts_top_20_features:both:5", "featcount_dist:both:5")
out <- table(unlist(lapply(strsplit(vars,":"), function(f){which(isOutlier(sce[[f[1]]], log=FALSE, nmads=as.numeric(f[3]), type=f[2]))})))
out <- as.numeric(names(out)[which(out > 0)])
fc <- data.frame(Feature=colData(sce)$log10_total_features, Count=colData(sce)$log10_total_counts, Group="Pass")
rownames(fc) <- rownames(colData(sce))
if (length(out) > 0) fc$Group[out] <- "Filter"
fc <- fc[order(fc$Group),]
if (length(out) > 0) sce <- sce[, -out]
sce <- as.Seurat(sce)
names(sce@assays) <- "RNA"
DefaultAssay(sce) <- "RNA"
pd <- ggplot(fc, aes(x=Feature, y=Count, colour=Group))+geom_point(stroke=0, shape=16, size=2, alpha=0.5)+
	stat_smooth(method=loess, se=FALSE, colour="black")+
	labs(title=NULL, x="Feature (log10)", y="Count (log10)")+
	scale_colour_manual(values=col_list[c(3, 1)])+
	guides(colour=guide_legend(override.aes=list(size=5)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, C=pd, nrow=1, widths=c(1, 1, 1, 7))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", nrow(mt_info)-length(which(dbi_info != "singlet")), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=12, height=6, dpi=200, "sce_qc_1.png")
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, C=pd, nrow=1, widths=c(1, 1, 1, 7))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", nrow(mt_info)-length(which(dbi_info != "singlet")), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=12, height=6, dpi=200, "sce_qc_1.pdf")

sce <- SCTransform(sce)
#ggsave(plot=LabelPoints(plot=VariableFeaturePlot(sce), points=head(VariableFeatures(sce), 10), repel=T), 
#	width=6, height=4, dpi=200, "rna_qc2.png")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:20)
sce <- FindNeighbors(sce, reduction="umap", dims=1:2)
sce <- FindClusters(sce, algorithm=4, method="igraph", resolution=0.2)
#sce <- FindNeighbors(sce)
#sce <- FindClusters(sce, algorithm=4, method="igraph")
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=5)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=wrap_plots(list(pa, pb), nrow=1), width=12, height=5, dpi=200, "rna_umap.png")

types <- names(table(sce[["seurat_clusters"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce, ident.1=x, group.by="seurat_clusters")})
names(types) <- names(table(sce[["seurat_clusters"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "rna_marker_cls.rds")

marker_info <- data.frame()
for (i in 1:length(types))
{
	terms <- types[[i]]
	terms$cls <- i
	marker_info <- rbind(marker_info, terms[1:min(nrow(terms), 10),])
}
write.csv(marker_info, "marker_info_cls.csv", quote=F)

#sce <- readRDS("patient.rds")
#types <- readRDS("rna_marker_cls.rds")
cell_marker <- read.csv("Cell_marker_Seq.csv", h=T)
cell_marker <- cell_marker[which(cell_marker$species == "Human" & cell_marker$cell_type == "Normal cell" & 
	cell_marker$tissue_type == "Blood" & cell_marker$Symbol != ""),]
genes <- intersect(rownames(sce), cell_marker$Symbol)
cell_marker$Pass <- "F"
for (gene in genes) cell_marker$Pass[which(cell_marker$Symbol == gene)] <- "T"
cell_marker <- cell_marker[which(cell_marker$Pass == "T"),]
ref_info <- table(cell_marker$cell_name)
ref_info <- tibble(type=names(ref_info[which(ref_info > 10)]), gene=list(""))
for (i in 1:nrow(ref_info)) ref_info$gene[i] <- list(unique(cell_marker$Symbol[which(cell_marker$cell_name == ref_info$type[i])]))

pred <- data.frame(Type="", FDR=0, Marker="", Other="", FDR.other="")
for (terms in types)
{
	term_list <- rownames(terms)
	find <- FALSE
	for (num in min(20, length(term_list)):length(term_list))
	{
		t <- enricher(term_list[1:num], TERM2GENE=ref_info, minGSSize=1)
		if (length(t$ID) > 0)
		{
			find <- TRUE
			pred <- rbind(pred, c(t$ID[1], t$p.adjust[1], gsub("/", ", ", t$geneID[1]), paste(t$ID[-1], collapse=", "), paste(t$p.adjust[-1], collapse=", ")))
			break
		}
	}
	if (find == FALSE) pred <- rbind(pred, c("NA", 1, "NA", "NA", 1))
}
pred <- pred[-1,]
rownames(pred) <- 1:nrow(pred)

sce$cell.type <- ""
for (i in 1:length(types)) sce$cell.type[which(sce$seurat_clusters == i)] <- pred$Type[i]
#sce$cell.type[which(sce$seurat_clusters == 12)] <- "Neuron"
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=4)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=4)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
pc <- UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=4)+
	labs(title="Type", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap.png")

saveRDS(sce, "patient.rds")

#sce <- readRDS("patient.rds")
terms <- c("CD19", "MS4A1", "HLA-DRA", "HLA-DRB1", "IGHM", "IGHD", "FCER2", 
	"CXCR4", "CD24", "CD86", "CD27", "CD38", "TNFRSF17", "BANK1", "MZB1", "MS4A1", "CD79A", 
	"TCL1A", "SOX4", "CD34", "CD3D", "CD3E", "CD14", "S100A9", "CD27", "CD83")
terms <- intersect(terms, rownames(sce))
terms <- c(terms, "nFeature_RNA", "nCount_RNA", "cell.mt")
pls <- lapply(terms, function(x)
{
	p <- FeaturePlot(sce, features=x, slot="counts", min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=x, x="UMAP1", y="UMAP2", colour=NULL)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pls, ncol=4), width=18, height=ceiling(length(terms)/4)*4, dpi=200, filename="rna_feature.png", limitsize=F)

terms <- c("BANK1", "MZB1", "MS4A1", "CD79A", "TCL1A", "SOX4", "CD34", "CD3D", "CD3E", "CD14", "S100A9", "CD27", "CD83")
terms <- intersect(terms, rownames(sce))
ggsave(plot=DotPlot(sce, features=terms, group.by="seurat_clusters", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=6, height=6, dpi=200, "cell_marker.png")

terms <- rev(c("CD19", "MS4A1", "HLA-DRA", "HLA-DRB1", "IGHM", "IGHD", "FCER2", 
	"CXCR4", "CD24", "CD86", "CD27", "CD38", "TNFRSF17"))
terms <- intersect(terms, rownames(sce))
ggsave(plot=DotPlot(sce, features=terms, group.by="seurat_clusters", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=6, height=6, dpi=200, "cell_marker2.png")

sce[["cell.type"]] <- "BMemory"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 2), 1] <- "BNaive"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 3), 1] <- "BNaive"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 5), 1] <- "BNaive"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 11), 1] <- "BPlasma"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 10), 1] <- "Mono"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 12), 1] <- "TC"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 14), 1] <- "MKC"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 15), 1] <- "DC"
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=4)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=4)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
pc <- UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=4)+
	labs(title="Type", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap.png")
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap.pdf")

types <- names(table(sce[["cell.type"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce, ident.1=x, group.by="cell.type")})
names(types) <- names(table(sce[["cell.type"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "rna_marker_type.rds")

marker_info <- data.frame()
for (i in names(types))
{
	terms <- types[[i]]
	terms$type <- i
	marker_info <- rbind(marker_info, terms[1:min(nrow(terms), 10),])
}
write.csv(marker_info, "marker_info_type.csv", quote=F)

saveRDS(sce, "patient.rds")

celltype <- c("BNaive", "BMemory", "BPlasma", "Mono", "DC", "MKC", "TC")
cell_marker <- read.delim("cell_marker.txt")
sce$cell.type <- factor(sce$cell.type, levels=celltype)
ggsave(plot=DotPlot(sce, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=5, height=6, dpi=200, "cell_marker_type.png")
ggsave(plot=DotPlot(sce, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=5, height=6, dpi=200, "cell_marker_type.pdf")


terms <- cell_marker$marker
terms <- intersect(terms, rownames(sce))
terms <- c(terms, "nFeature_RNA", "nCount_RNA", "cell.mt")
pls <- lapply(terms, function(x)
{
	t <- x
	if (t != "nFeature_RNA" & t != "nCount_RNA" & t != "cell.mt") t <- paste0(x, "(", cell_marker$type[match(x, cell_marker$marker)], ")")
	p <- FeaturePlot(sce, features=x, slot="counts", min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=t, x="UMAP1", y="UMAP2", colour=NULL)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature.png", limitsize=F)
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature.pdf", limitsize=F)



terms <- cell_marker$marker
terms <- intersect(terms, rownames(sce))
terms <- c(terms, "nFeature_RNA", "nCount_RNA", "cell.mt")
pls <- lapply(terms, function(x)
{
	t <- x
	if (t != "nFeature_RNA" & t != "nCount_RNA" & t != "cell.mt") t <- paste0(x, "(", cell_marker$type[match(x, cell_marker$marker)], ")")
	p <- FeaturePlot(sce, features=x, slot="counts", min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=t, x="UMAP1", y="UMAP2", colour=NULL)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature.png", limitsize=F)
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature.pdf", limitsize=F)

##############################################
#BCR-seq--Normal
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
library(ggbeeswarm)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", 
	"#9e4747", "#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8")
text_size <- 11
title_size <- 13
tag_thm <- theme(plot.tag=element_text(size=14, colour="black", face="bold"), plot.margin=margin())


#fileList <- list.files(path="data/normal", pattern="*.counts.txt$")
#features <- rownames(read.delim(paste0("data/normal/", fileList[1])))
#for(file in fileList[-1])
#{
#	data_raw <- read.delim(paste0("data/normal/", file))
#	features <- intersect(features, rownames(data_raw))
#}
#mat_norm <- read.delim(paste0("data/normal/", fileList[1]))[features,]
#for(file in fileList[-1])
#{
#	data_raw <- read.delim(paste0("data/normal/", file))[features,]
#	mat_norm <- cbind(mat_norm, data_raw)
#}
#write.table(mat_norm, "normal_rna.txt", sep="\t", )
#fileList <- list.files(path="data/normal", pattern="*_annotations.csv$")
#mat_bcr <- data.frame()
#for(file in fileList)
#{
#	data_raw <- read.csv(paste0("data/normal/", file), h=T)
#	data_raw$barcode <- paste0(gsub("_.*", ".", gsub(".*_HC", "HC", file)), gsub("-1", "", data_raw$barcode))
#	mat_bcr <- rbind(mat_bcr, data_raw)
#}
#write.csv(mat_bcr, "data/normal_bcr.csv", row.names=F, quote=F)
sce <- CreateSeuratObject(read.delim("data/normal_rna.txt"), project="normal")
sce[["cell.sample"]] <- Idents(sce)
sce[["cell.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
cs_info <- data.frame(sce@meta.data)
pa <- ggplot(cs_info, aes(x="Term", y=log10(nFeature_RNA)))+geom_violin(trim=F, fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+scale_y_continuous(breaks=seq(2.8, 3.8, 0.2))+
	labs(title=NULL, x="Genes (log10)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pb <- ggplot(cs_info, aes(x="Term", y=log10(nCount_RNA)))+geom_violin(trim=F, fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Counts (log10)", y=NULL)+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pc <- ggplot(cs_info, aes(x="Term", y=cell.mt))+geom_violin(trim=F, fill=col_list[1], color=col_list[1])+
	geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_y_continuous(breaks=seq(5, 25, 5), limit=c(0, 30))+scale_x_discrete(breaks=NULL)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
ggsave(plot=pa+pb+pc, width=6, height=4, dpi=200, "rna_qc_norm.png")
ggsave(plot=pa+pb+pc, width=6, height=4, dpi=200, "rna_qc_norm.pdf")


sce <- as.SingleCellExperiment(sce)
id <- grep("^IG[LHK].*", rownames(sce))
if (length(id) > 0) sce <- sce[-id,]
#sce <- scDblFinder(sce, samples="cell.sample")
sce$scDblFinder.class <- "singlet"
dbi_info <- sce$scDblFinder.class
mt_limit <- 15
mt_info <- data.frame(value=sce$cell.mt, type="F")
mt_info$type[which(mt_info$value < mt_limit)] <- "P"
cs_info <- data.frame(value=colSums(counts(sce) > 0), umi=colSums(counts(sce)), type="F")
cs_limit <- c(800, 3000)
cs_info$type[which(cs_info$value > cs_limit[1] & cs_info$umi > cs_limit[1]*2 & cs_info$value < cs_limit[2] & cs_info$umi < cs_limit[2]*3)] <- "P"
sce <- sce[, which(mt_info$type == "P" & cs_info$type == "P" & sce$scDblFinder.class == "singlet")]
sce <- sce[rowSums(counts(sce) > 0) > 10,]
pa <- ggplot(mt_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Mito. (%)", y=NULL)+scale_y_continuous(breaks=seq(10, 90, 10))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=mt_limit)+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pb <- ggplot(cs_info, aes(x="Term", y=umi, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Counts", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1]*2, seq(5000, 50000, 5000)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=c(cs_limit[1]*2, cs_limit[2]*3))+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())
pc <- ggplot(cs_info, aes(x="Term", y=value, colour=type))+geom_quasirandom(size=0.5, groupOnX=T)+
	scale_colour_manual(values=col_list[c(3, 1)])+guides(colour="none")+
	labs(title=NULL, x="Feature", y=NULL)+scale_y_continuous(breaks=c(cs_limit[1], seq(500, 5000, 500)))+
	scale_x_discrete(breaks=NULL)+geom_hline(yintercept=c(cs_limit[1], cs_limit[2]))+
	theme(axis.text=element_text(colour="black"), axis.line=element_line(colour="black"), panel.background=element_blank())

sce <- addPerCellQC(sce, percent_top=c(20, 50, 100, 200))
sce$total_features <- sce$detected
sce$log10_total_features <- log10(sce$detected)
sce$total_counts <- sce$sum
sce$log10_total_counts <- log10(sce$sum)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features
mod <- loess(colData(sce)$log10_total_features~colData(sce)$log10_total_counts)
pred <- predict(mod, newdata=data.frame(log10_total_counts=colData(sce)$log10_total_counts))
sce$featcount_dist <- colData(sce)$log10_total_features - pred
sce$pct_counts_top_20_features <- colData(sce)[[intersect(c("percent_top_20","pct_counts_in_top_20_features","percent.top_20"), colnames(colData(sce)))[[1]]]]
sce$pct_counts_top_50_features <- colData(sce)[[intersect(c("percent_top_50","pct_counts_in_top_50_features","percent.top_50"), colnames(colData(sce)))[[1]]]]
vars <- c("log10_total_counts:both:5", "log10_total_features:both:5", "pct_counts_top_20_features:both:5", "featcount_dist:both:5")
out <- table(unlist(lapply(strsplit(vars,":"), function(f){which(isOutlier(sce[[f[1]]], log=FALSE, nmads=as.numeric(f[3]), type=f[2]))})))
out <- as.numeric(names(out)[which(out > 0)])
fc <- data.frame(Feature=colData(sce)$log10_total_features, Count=colData(sce)$log10_total_counts, Group="Pass")
rownames(fc) <- rownames(colData(sce))
if (length(out) > 0) fc$Group[out] <- "Filter"
fc <- fc[order(fc$Group),]
if (length(out) > 0) sce <- sce[, -out]
sce <- as.Seurat(sce)
names(sce@assays) <- "RNA"
DefaultAssay(sce) <- "RNA"
pd <- ggplot(fc, aes(x=Feature, y=Count, colour=Group))+geom_point(stroke=0, shape=16, size=2, alpha=0.5)+
	stat_smooth(method=loess, se=FALSE, colour="black")+
	labs(title=NULL, x="Feature (log10)", y="Count (log10)")+
	scale_colour_manual(values=col_list[c(3, 1)])+
	guides(colour=guide_legend(override.aes=list(size=5)))+
	theme(axis.line=element_line(linetype=1, colour='black'), 
	legend.key=element_blank(), legend.background=element_blank(), 
	panel.background=element_rect(0, linetype=0))
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, C=pd, nrow=1, widths=c(1, 1, 1, 7))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", nrow(mt_info)-length(which(dbi_info != "singlet")), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=12, height=6, dpi=200, "sce_qc_1_norm.png")
ggsave(plot=wrap_plots(A=pa, B=pb, C=pc, C=pd, nrow=1, widths=c(1, 1, 1, 7))+
	plot_annotation(title=paste0("The QC results of scRNA-seq (", dim(sce)[2], "/", nrow(mt_info)-length(which(dbi_info != "singlet")), ")"), 
	theme=theme(plot.title=element_text(size=16, hjust=0.5))), 
	width=12, height=6, dpi=200, "sce_qc_1_norm.pdf")

sce <- SCTransform(sce)
#ggsave(plot=LabelPoints(plot=VariableFeaturePlot(sce), points=head(VariableFeatures(sce), 10), repel=T), 
#	width=6, height=4, dpi=200, "rna_qc2.png")
sce <- RunPCA(sce)
sce <- RunUMAP(sce, dims=1:5)
sce <- FindNeighbors(sce, reduction="umap", dims=1:2)
sce <- FindClusters(sce, algorithm=4, method="igraph", resolution=0.2)
#sce <- FindNeighbors(sce)
#sce <- FindClusters(sce, algorithm=4, method="igraph")
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=5)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=5)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=wrap_plots(list(pa, pb), nrow=1), width=12, height=5, dpi=200, "rna_umap_norm.png")

types <- names(table(sce[["seurat_clusters"]][,1]))
types <- lapply(types, function(x){FindMarkers(sce, ident.1=x, group.by="seurat_clusters")})
names(types) <- names(table(sce[["seurat_clusters"]][,1]))
for (i in 1:length(types)) types[[i]] <- types[[i]][which(types[[i]]$p_val_adj < 0.01 & types[[i]]$avg_log2FC > 0.1),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$avg_log2FC, decreasing=T),]
for (i in 1:length(types)) types[[i]] <- types[[i]][order(types[[i]]$p_val_adj),]
saveRDS(types, "rna_marker_cls_norm.rds")

marker_info <- data.frame()
for (i in 1:length(types))
{
	terms <- types[[i]]
	terms$cls <- i
	marker_info <- rbind(marker_info, terms[1:min(nrow(terms), 10),])
}
write.csv(marker_info, "marker_info_cls_norm.csv", quote=F)

#types <- readRDS("rna_marker_cls_norm.rds")
cell_marker <- read.csv("Cell_marker_Seq.csv", h=T)
cell_marker <- cell_marker[which(cell_marker$species == "Human" & cell_marker$cell_type == "Normal cell" & 
	cell_marker$tissue_type == "Blood" & cell_marker$Symbol != ""),]
genes <- intersect(rownames(sce), cell_marker$Symbol)
cell_marker$Pass <- "F"
for (gene in genes) cell_marker$Pass[which(cell_marker$Symbol == gene)] <- "T"
cell_marker <- cell_marker[which(cell_marker$Pass == "T"),]
ref_info <- table(cell_marker$cell_name)
ref_info <- tibble(type=names(ref_info[which(ref_info > 10)]), gene=list(""))
for (i in 1:nrow(ref_info)) ref_info$gene[i] <- list(unique(cell_marker$Symbol[which(cell_marker$cell_name == ref_info$type[i])]))

pred <- data.frame(Type="", FDR=0, Marker="", Other="", FDR.other="")
for (terms in types)
{
	term_list <- rownames(terms)
	find <- FALSE
	for (num in min(20, length(term_list)):length(term_list))
	{
		t <- enricher(term_list[1:num], TERM2GENE=ref_info, minGSSize=1)
		if (length(t$ID) > 0)
		{
			find <- TRUE
			pred <- rbind(pred, c(t$ID[1], t$p.adjust[1], gsub("/", ", ", t$geneID[1]), paste(t$ID[-1], collapse=", "), paste(t$p.adjust[-1], collapse=", ")))
			break
		}
	}
	if (find == FALSE) pred <- rbind(pred, c("NA", 1, "NA", "NA", 1))
}
pred <- pred[-1,]
rownames(pred) <- 1:nrow(pred)

sce$cell.type <- ""
for (i in 1:length(types)) sce$cell.type[which(sce$seurat_clusters == i)] <- pred$Type[i]
#sce$cell.type[which(sce$seurat_clusters == 12)] <- "Neuron"
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
pc <- UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Type", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap_norm.png")

#saveRDS(sce, "normal.rds")

#sce <- readRDS("normal.rds")
terms <- c("CD19", "MS4A1", "HLA-DRA", "HLA-DRB1", "IGHM", "IGHD", "FCER2", 
	"CXCR4", "CD24", "CD86", "CD27", "CD38", "TNFRSF17", "BANK1", "MZB1", "MS4A1", "CD79A", 
	"TCL1A", "SOX4", "CD34", "CD3D", "CD3E", "CD14", "S100A9", "CD27", "CD83")
terms <- intersect(terms, rownames(sce))
terms <- c(terms, "nFeature_RNA", "nCount_RNA", "cell.mt")
pls <- lapply(terms, function(x)
{
	p <- FeaturePlot(sce, features=x, slot="counts", min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=x, x="UMAP1", y="UMAP2", colour=NULL)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pls, ncol=4), width=18, height=ceiling(length(terms)/4)*4, dpi=200, filename="rna_feature_norm.png", limitsize=F)

terms <- c("BANK1", "MZB1", "MS4A1", "CD79A", "TCL1A", "SOX4", "CD34", "CD3D", "CD3E", "CD14", "S100A9", "CD27", "CD83")
terms <- intersect(terms, rownames(sce))
ggsave(plot=DotPlot(sce, features=terms, group.by="seurat_clusters", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=6, height=6, dpi=200, "cell_marker_norm.png")

terms <- rev(c("CD19", "MS4A1", "HLA-DRA", "HLA-DRB1", "IGHM", "IGHD", "FCER2", 
	"CXCR4", "CD24", "CD86", "CD27", "CD38", "TNFRSF17"))
terms <- intersect(terms, rownames(sce))
ggsave(plot=DotPlot(sce, features=terms, group.by="seurat_clusters", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=6, height=6, dpi=200, "cell_marker2_norm.png")

sce[["cell.type"]] <- "BNaive"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 2), 1] <- "BMemory"
sce[["cell.type"]][which(sce[["seurat_clusters"]][,1] == 6), 1] <- "BMemory"
pa <- UMAPPlot(sce, group.by="cell.sample", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Sample", x="UMAP1", y="UMAP2", colour=NULL)
pb <- UMAPPlot(sce, group.by="seurat_clusters", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Cluster", x="UMAP1", y="UMAP2", colour=NULL)
pc <- UMAPPlot(sce, group.by="cell.type", pt.size=1, label=T, label.size=4, cols=col_list)+
	labs(title="Type", x="UMAP1", y="UMAP2", colour=NULL)
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap_norm.png")
ggsave(plot=pa+pb+pc, width=28, height=8, dpi=200, "rna_umap_norm.pdf")

saveRDS(sce, "normal.rds")

celltype <- c("BNaive", "BMemory", "BPlasma", "Mono", "DC", "MKC", "TC")
cell_marker <- read.delim("cell_marker.txt")
sce$cell.type <- factor(sce$cell.type, levels=celltype)
ggsave(plot=DotPlot(sce, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=3, height=6, dpi=200, "cell_marker_type_norm.png")
ggsave(plot=DotPlot(sce, features=rev(cell_marker$marker), group.by="cell.type", 
	dot.min=0.1, col.min=0.1, cols=c(brewer.pal(11,"Spectral")[6], brewer.pal(11,"Spectral")[11]))+coord_flip()+
	labs(title=NULL, x=NULL, y=NULL)+
	guides(color=guide_colorbar(title='Avg Expr'), size=guide_legend(title="Per Expr"))+
	theme(axis.text.x=element_text(size=text_size, angle=270, vjust=0.5, hjust=0, colour="black"), 
	axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
	axis.text.y=element_text(size=text_size, colour="black", vjust=0.5, hjust=1), 
	legend.title=element_text(size=title_size, colour="black", face="bold"), 
	legend.text=element_text(size=text_size, colour="black")), 
	width=3, height=6, dpi=200, "cell_marker_type_norm.pdf")

terms <- cell_marker$marker
terms <- intersect(terms, rownames(sce))
terms <- c(terms, "nFeature_RNA", "nCount_RNA", "cell.mt")
pls <- lapply(terms, function(x)
{
	t <- x
	if (t != "nFeature_RNA" & t != "nCount_RNA" & t != "cell.mt") t <- paste0(x, "(", cell_marker$type[match(x, cell_marker$marker)], ")")
	p <- FeaturePlot(sce, features=x, slot="counts", min.cutoff=0, cols=c("gray98", "red2"))+
		labs(title=t, x="UMAP1", y="UMAP2", colour=NULL)+
		theme(plot.title=element_text(size=14, hjust=0.5, face="bold"), panel.background=element_blank(), 
		axis.line=element_line(colour="black"), axis.ticks=element_blank(), axis.text=element_blank(), 
		axis.title=element_text(colour="black", size=14), legend.title=element_text(colour="black", size=14), 
		legend.text=element_text(colour="black", size=12))
	p$data <- p$data[order(p$data[, 4]),]
	return(p)
})
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature_norm.png", limitsize=F)
ggsave(plot=wrap_plots(pls, ncol=5), width=22, height=ceiling(length(terms)/5)*4, dpi=200, filename="rna_feature_norm.pdf", limitsize=F)












#################################################################################
Fig3A-F
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
library(ggbeeswarm)
options(stringsAsFactors=FALSE)

col_list <- c("#46998b", "#847acc", "#ef8560", "#6994b3", "#d1934b", "#8fb350", "#de9cba", "#7b469e", "#9e4747", 
	"#1e8751", "#cc9a04", "#4bb35b", "#e13344", "#855949", "#3b4992", "#6e84b8", "#9AD3DD", "#3071A9", "#31B3C2")
text_size <- 11
title_size <- 13
tag_thm <- theme(plot.tag=element_text(size=20, colour="black"), plot.margin=margin())

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



sce <- sce_patient
bcrinfo <- bcrinfo_patient
bcrinfo <- bcrinfo[which(bcrinfo$c_gene != "" & bcrinfo$c_gene != "IGKC" & bcrinfo$c_gene != "None"),]
if (length(grep("IGLC", bcrinfo$c_gene)) > 0) bcrinfo <- bcrinfo[-grep("IGLC", bcrinfo$c_gene),]
celltype <- c("BNaive", "BMemory", "BPlasma")
bcs <- c()
terms <- c()
for (type in celltype) bcs <- c(bcs, colnames(sce)[which(sce$cell.type == type)])
bcs <- intersect(bcs, bcrinfo$barcode)
sce <- subset(sce, cells=bcs)
sce[["cell.cg"]] <- gsub("[0-9]+$", "", bcrinfo$c_gene[match(colnames(sce), bcrinfo$barcode)])
sce[["cell.clone"]] <- "Non-clonal"
sce$cell.clone[match(intersect(bcrinfo$barcode[which(bcrinfo$cc > 1)], colnames(sce)), colnames(sce))] <- "Clonal"
sce$cell.clone <- factor(sce$cell.clone, levels=c("Clonal", "Non-clonal"))
sce_a <- sce
paa <- UMAPPlot(sce_a, group.by="cell.type", pt.size=0.5, label=F, order=c("BNaive", "BMemory", "BPlasma"), cols=col_list[c(19,18,17)])+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour="Type")+theme(legend.position="none", plot.margin=margin())
pba <- UMAPPlot(sce_a, group.by="cell.cg", pt.size=0.5, label=F, order=c("IGHA", "IGHG", "IGHM", "IGHD"), cols=col_list[c(4:1)])+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour="CGene")+theme(legend.position="none", plot.margin=margin())
pca <- UMAPPlot(sce_a, group.by="cell.clone", pt.size=0.5, label=F, order=c("Clonal", "Non-clonal"), cols=col_list[c(4, 7)])+
	labs(title=NULL, x="UMAP1", y="UMAP2", colour="Clonal")+theme(legend.position="none", plot.margin=margin())

sce <- sce_normal
cell_normal <- data.frame(table(sce$cell.type), Group="Normal")
sce <- sce_patient
cell_patient <- data.frame(table(paste0(sce$cell.sample, "+", sce$cell.type)), Group="Patient")
cell_patient$Group <- gsub("\\+.*$", "", cell_patient$Var1)
cell_patient$Var1 <- gsub("^.*\\+", "", cell_patient$Var1)
res <- rbind(cell_normal, cell_patient)
ids <- c()
for (type in c("BNaive", "BMemory", "BPlasma")) ids <- c(ids, which(res$Var1 == type))
res <- res[ids,]
res$Rate <- res$Freq
for (term in unique(res$Group))
{
	ids <- which(res$Group ==  term)
	res$Rate[ids] <- res$Rate[ids]*100/sum(res$Freq[ids])
}
res$Var1 <- factor(res$Var1, levels=rev(c("BNaive", "BMemory", "BPlasma")))
resa <- res
pab <- ggplot(resa, aes(x=Group, y=Rate, fill=Var1))+geom_bar(stat="identity")+
	scale_y_continuous(breaks=seq(0,100,25), labels=paste0(seq(0,100,25), "%"), expand=c(0, 0))+
	scale_fill_manual(values=col_list[c(19,18,17)])+labs(title=NULL, fill="Type", x=NULL, y=NULL)+
	theme(panel.background=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, angle=-45, hjust=0, vjust=0.5, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks.x=element_blank(), axis.line=element_line(colour="black"), plot.margin=margin())

sce <- sce_normal
bcrinfo <- bcrinfo_normal
bcrinfo <- bcrinfo[which(bcrinfo$c_gene != "" & bcrinfo$c_gene != "IGKC" & bcrinfo$c_gene != "None"),]
if (length(grep("IGLC", bcrinfo$c_gene)) > 0) bcrinfo <- bcrinfo[-grep("IGLC", bcrinfo$c_gene),]
celltype <- c("BNaive", "BMemory", "BPlasma")
res_total <- data.frame()
for (type in celltype)
{
	ids <- c()
	for (bc in colnames(sce)[which(sce$cell.type == type)]) ids <- c(ids, which(bcrinfo$barcode == bc))
	if (length(ids) == 0) next
	res <- data.frame(table(gsub("[0-9]+$", "", bcrinfo$c_gene[ids])))
	res$Rate <- res$Freq*100/sum(res$Freq)
	res$Type <- type
	res_total <- rbind(res_total, res)
}
res_total$Type <- factor(res_total$Type, levels=celltype)
res_total$Sample <- "Normal"
res_normal <- res_total
sce <- sce_patient
sce$cell.temp <- paste0(sce$cell.sample, "+", sce$cell.type)
bcrinfo <- bcrinfo_patient
bcrinfo <- bcrinfo[which(bcrinfo$c_gene != "" & bcrinfo$c_gene != "IGKC" & bcrinfo$c_gene != "None"),]
if (length(grep("IGLC", bcrinfo$c_gene)) > 0) bcrinfo <- bcrinfo[-grep("IGLC", bcrinfo$c_gene),]
celltype <- c("BNaive", "BMemory", "BPlasma")
res_total <- data.frame()
terms <- unique(paste0(sce$cell.sample, "+", celltype))
for (type in terms)
{
	ids <- c()
	for (bc in colnames(sce)[which(sce$cell.temp == type)]) ids <- c(ids, which(bcrinfo$barcode == bc))
	if (length(ids) == 0) next
	res <- data.frame(table(gsub("[0-9]+$", "", bcrinfo$c_gene[ids])))
	res$Rate <- res$Freq*100/sum(res$Freq)
	res$Type <- type
	res_total <- rbind(res_total, res)
}
res_total$Sample <- gsub("\\+.*$", "", res_total$Type)
res_total$Type <- factor(gsub("^.*\\+", "", res_total$Type), levels=celltype)
res_patient <- res_total
res_total <- rbind(res_normal, res_patient)
terms <- as.character(unique(res_total$Var1))
terms <- terms[order(terms)]
terms <- c("IGHA", "IGHG", "IGHM", "IGHD")
res_total$Var1 <- factor(res_total$Var1, levels=terms)
resb <- res_total
pbb <- ggplot(resb, aes(x=Type, y=Rate, fill=Var1))+geom_bar(stat="identity")+
	scale_y_continuous(breaks=seq(0,100,25), labels=paste0(seq(0,100,25), "%"), expand=c(0, 0))+
	scale_fill_manual(values=rev(col_list[c(4:1)]))+
	labs(title=NULL, fill="Type", x=NULL, y=NULL)+facet_wrap(. ~ Sample, nrow=1)+
	theme(panel.background=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, angle=-45, hjust=0, vjust=0.5, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(fill="white", colour="white"), 
	axis.ticks.x=element_blank(), axis.line=element_line(colour="black"), plot.margin=margin())



###########################################################################3
Fig3G
cs_info <- data.frame(Cell=colnames(sce), Clone=sce$cell.clone, Sample=sce$cell.sample, Type=sce$cell.type, Value=0)
cs_info <- cs_info[match(intersect(cs_info$Cell, colnames(sce_add)), cs_info$Cell),]

genes <- rownames(sce_add)[grep("^IGHV3-", rownames(sce_add))]
genes <- genes[which(rowSums(sce_add[["RNA"]]@counts[genes,]) > 10)]
#sce_add <- AddModuleScore(sce_add, features=genes, name="IGHV", nbin=23)
sce_add$IGHV1 <- colSums(sce_add[["RNA"]]@counts[genes,])
cs_info$Value <- sce_add$IGHV1[match(cs_info$Cell, colnames(sce_add))]
cs_info$Type <- factor(cs_info$Type, levels=c("BNaive", "BMemory", "BPlasma"))
cs_info$Clone <- factor(cs_info$Clone, levels=c("Non-clonal", "Clonal"))
cs_info_c <- cs_info
pcc <- ggplot(cs_info_c, aes(x=Type, y=Value, fill=Clone, colour=Clone))+geom_violin(position=position_dodge(0.9))+
	geom_boxplot(width=0.6, outlier.alpha=0, fill="white", position=position_dodge(0.9))+facet_wrap(~Sample, nrow=1)+
	scale_fill_manual(values=col_list[c(4, 7)], drop=F)+
	scale_colour_manual(values=col_list[c(4, 7)], drop=F)+
	#geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x=NULL, y="Expression of IGHV")+
	theme(panel.background=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, angle=-45, hjust=0, vjust=0.5, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(fill="white", colour="white"), 
	axis.ticks.x=element_blank(), axis.line=element_line(colour="black"), plot.margin=margin())

freq_mtx <- read.delim("freq_mtx_ighv3.tsv")
#freq_mtx <- freq_mtx[grep("-A-G$", rownames(freq_mtx)),]
cs <- colSums(freq_mtx > 0)
names(cs) <- gsub("\\.", "-", names(cs))
cs_info <- data.frame(Cell=colnames(sce), Clone=sce$cell.clone, Sample=sce$cell.sample, Type=sce$cell.type, Value=0)
cs_info <- cs_info[match(intersect(cs_info$Cell, names(cs)), cs_info$Cell),]
cs_info$Value <- cs[match(cs_info$Cell, names(cs))]
cs_info$Type <- factor(cs_info$Type, levels=c("BNaive", "BMemory", "BPlasma"))
cs_info$Clone <- factor(cs_info$Clone, levels=c("Non-clonal", "Clonal"))
cs_info_d <- cs_info
pcd <- ggplot(cs_info_d, aes(x=Type, y=Value, fill=Clone, colour=Clone))+geom_violin(position=position_dodge(0.9))+
	geom_boxplot(width=0.2, outlier.alpha=0, fill="white", position=position_dodge(0.9))+facet_wrap(~Sample, nrow=1)+
	scale_fill_manual(values=col_list[c(4, 7)], drop=F)+
	scale_colour_manual(values=col_list[c(4, 7)], drop=F)+
	#geom_boxplot(fill="white", outlier.alpha=0, width=0.3)+
	labs(title=NULL, x=NULL, y="Mutations of IGHV")+
	theme(panel.background=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	axis.text.x=element_text(size=text_size, angle=-45, hjust=0, vjust=0.5, colour="black"), 
	plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(fill="white", colour="white"), 
	axis.ticks.x=element_blank(), axis.line=element_line(colour="black"), plot.margin=margin())

ggsave(plot=wrap_plots(list(wrap_elements(paa)+tag_thm, wrap_elements(pab)+tag_thm, wrap_elements(pba)+tag_thm, 
	wrap_elements(pbb)+tag_thm), nrow=1, widths=c(5,5,5,10))+plot_annotation(tag_levels="A"), 
	width=18, height=4, dpi=200, filename="bcr_fig3_a.png", limitsize=F)
ggsave(plot=wrap_plots(list(wrap_elements(paa)+tag_thm, wrap_elements(pab)+tag_thm, wrap_elements(pba)+tag_thm, 
	wrap_elements(pbb)+tag_thm), nrow=1, widths=c(5,5,5,10))+plot_annotation(tag_levels="A"), 
	width=18, height=4, dpi=200, filename="bcr_fig3_a.pdf", limitsize=F)
ggsave(plot=wrap_plots(list(wrap_elements(pca+tag_thm+plot_annotation(tag_levels=list(c("E")))), 
	wrap_elements(wrap_plots(list(pcb, pcc, pcd), widths=c(1,4,4))+plot_layout(guides="collect")+
	tag_thm+plot_annotation(tag_levels=list(c("F", "G", "H"))))), nrow=1, widths=c(5,20))+plot_layout(guides="collect"), 
	width=18, height=4, dpi=200, filename="bcr_fig3_b.png", limitsize=F)
ggsave(plot=wrap_plots(list(wrap_elements(pca+tag_thm+plot_annotation(tag_levels=list(c("E")))), 
	wrap_elements(wrap_plots(list(pcb, pcc, pcd), widths=c(1,4,4))+plot_layout(guides="collect")+
	tag_thm+plot_annotation(tag_levels=list(c("F", "G", "H"))))), nrow=1, widths=c(5,20))+plot_layout(guides="collect"), 
	width=18, height=4, dpi=200, filename="bcr_fig3_b.pdf", limitsize=F)






#######################################################################
Fig3H

library(Seurat)
library(Signac)
library(ggrepel)
library(dplyr)
library(dittoSeq)
library(MySeuratWrappers)
library(ggpubr)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Fig3/")
BCR=readRDS("/md01/xuyc7/work/BCR/bcr_expr.rds")
Idents(BCR)="cell.type"
DefaultAssay(BCR)="RNA"
cell=c("BNaive","BMemory","BPlasma")
names(cell) <- levels(BCR)
BCR <- RenameIdents(BCR, cell)
test_sign=list(c("Clonal","Non-clonal"))
gene=c("IRF1","IRF2","IRF4","TCF12","STAT1","STAT2",
       "IFNAR1","IFNAR2","TYK2","JAK1","JAK2","IL10RB","IFNLR1","STAT1","STAT2","IRF9","IFNGR1","IFNGR2")
pdf("Fig3G_TF_gene_exp_mean_vlnplot_RNA.pdf",width=6,height=4)
for(m in 1:length(gene)){
p=VlnPlot(BCR, features = gene[m], 
        split.by = 'cell.clone', pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols =c("#AA6F9EFF","#437DBFFF")) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("Cell") + ylab("Gene Expression") + 
       stat_compare_means( method = "wilcox.test") +
       stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 8, angle = 0))+
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
      print(p)
}
dev.off()

RNA_data=GetAssayData(BCR)
RNA_data1=RNA_data[gene,]
RES=c()
for(i in 1:length(gene)){
  for(n in 1:length(cell)){
  exp_Pat=RNA_data1[gene[i],which(BCR$cell.type==cell[n]&BCR$cell.clone=="Clonal")]
  exp_Nor=RNA_data1[gene[i],which(BCR$cell.type==cell[n]&BCR$cell.clone=="Non-clonal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  res=cbind(gene[i],as.character(cell[n]),Test_p,FC)
  RES=rbind(RES,res)
}
}

RES=c()
for(i in 1:length(gene)){
  exp_Pat=RNA_data1[gene[i],which(BCR$cell.clone=="Clonal")]
  exp_Nor=RNA_data1[gene[i],which(BCR$cell.clone=="Non-clonal")]
  Test=wilcox.test(exp_Pat,exp_Nor)
  Test_p=Test$p.value
  FC=mean(exp_Pat)/mean(exp_Nor)
  res=cbind(gene[i],as.character(cell[n]),Test_p,FC)
  RES=rbind(RES,res)
}

Idents(BCR)="cell.clone"
test_sign=list(c("Clonal","Non-clonal"))
gene=c("IRF1","IRF2","IRF4","TCF12","STAT1","STAT2",
       "IFNAR1","IFNAR2","TYK2","JAK1","JAK2","IL10RB","IFNLR1","STAT1","STAT2","IRF9","IFNGR1","IFNGR2",
       "IGHM","IGHG1","IGHG2","IGHG3","IGHG4","IGHD","IGHA1","IGHA2")

pdf("Fig3H_TF_gene_exp_mean_vlnplot_RNA_merge.pdf",width=18,height=24)
#for(m in 1:length(gene)){
VlnPlot(BCR, features = gene, 
        split.by = 'cell.clone',
         pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols =c("#AA6F9EFF","#437DBFFF")) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("Cell") + ylab("Gene Expression") + 
       stat_compare_means( method = "wilcox.test") +
       stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 8, angle = 0))+
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
dev.off()

IGHV=rownames(RNA_data)[grep("IGHV",rownames(RNA_data))]
pdf("Fig3H_IGHV_exp_mean_vlnplot_RNA_merge.pdf",width=36,height=36)
#for(m in 1:length(gene)){
VlnPlot(BCR, features = IGHV[101:146], 
        split.by = 'cell.clone',
         pt.size = 0.0, 
        combine = T, #split.plot = T,
        cols =c("#AA6F9EFF","#437DBFFF")) + 
       theme_bw() + 
       theme(legend.title = element_blank()) +
       stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = 0.6,
                 position = position_dodge(width = .70)) +   
       xlab("Cell") + ylab("Gene Expression") + 
       stat_compare_means( method = "wilcox.test") +
       stat_compare_means(comparisons = test_sign,size = 5, label = "p.signif")+
       theme(text = element_text(size = 8, angle = 0))+
       guides(fill = guide_legend(override.aes = list(linetype = 0)),
              color = guide_legend(override.aes = list(linetype = 0))
      )
dev.off()



#######################################################333
Fig 3I-J

sce <- sce_patient
bcrinfo <- bcrinfo_patient
pls <- list()
sample_list <- paste0("Patient", 1:5)
for (sample in sample_list)
{
	bcrinfo_sub <- bcrinfo[which(bcrinfo$sample == sample),]
	bcrinfo_sub <- bcrinfo_sub[which(bcrinfo_sub$type == "BMemory" | bcrinfo_sub$type == "BPlasma"),]
	res <- data.frame()
	for (term in unique(bcrinfo_sub$cdr3)) res <- rbind(res, data.frame(Term=term, 
		BM=length(unique(bcrinfo_sub$barcode[which(bcrinfo_sub$cdr3 == term & bcrinfo_sub$type == "BMemory")])), 
		BP=length(unique(bcrinfo_sub$barcode[which(bcrinfo_sub$cdr3 == term & bcrinfo_sub$type == "BPlasma")]))))
	ids <- which(res$BM > 0 & res$BP > 0)
	if (length(ids) == 0) next
	res <- res[ids,]
	res$Total <- res$BM + res$BP
	res <- res[order(res$Total, decreasing=T),]
	ids <- which(res$Total > 5)
	if (length(ids) == 0) next
	res <- res[ids,]
	rec <- data.frame()
	for (i in 1:nrow(res)) rec <- rbind(rec, data.frame(Group=c("BPlasma", "BMemory"), Term=res$Term[i], Value=c(res$BM[i], res$BP[i])))
	rec$Term <- factor(rec$Term, levels=rev(res$Term))
	rec$Group <- factor(rec$Group, levels=c("BPlasma", "BMemory"))
	pls[[sample]] <- ggplot(rec, aes(x=Group, y=Term, size=Value))+geom_point(color=col_list[2])+
		labs(title=sample, x=NULL, y=NULL)+
		theme(plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
		axis.text=element_text(size=text_size, colour="black"), panel.background=element_blank(), 
		legend.title=element_text(colour="black", face="bold"), 
		axis.ticks=element_blank(), axis.line=element_line(colour="black"))
}
ggsave(plot=wrap_plots(pls, ncol=1, heights=c(2,10,1)), width=6, height=12, dpi=200, "bm_mp.png", limitsize=F)

sample_list <- paste0("Patient", 1:4)
sce <- sce_patient
bcrinfo <- bcrinfo_patient
bcrinfo_sub <- bcrinfo[which(bcrinfo$type == "BPlasma"),]
res <- data.frame()
cdr_info <- data.frame()
for (term in unique(bcrinfo_sub$cdr3))
{
	ids <- which(bcrinfo_sub$cdr3 == term)
	if (length(ids) < 3) next
	info <- data.frame(table(gsub("^.*_", "", unique(paste0(bcrinfo_sub$barcode[ids], "_", bcrinfo_sub$sample[ids])))), Term=term)
	cdr_info <- rbind(cdr_info, data.frame(Term=term, Count=sum(info$Freq), Group=info$Var1[which.max(info$Freq)]))
	res <- rbind(res, info)
}
cdr_info$Group <- factor(cdr_info$Group, levels=sample_list)
cdr_info <- cdr_info[order(cdr_info$Count, decreasing=T),]
cdr_info <- cdr_info[order(cdr_info$Group),]
cdr_info$SN <- 0
for (i in 1:nrow(cdr_info)) cdr_info$SN[i] <- length(unique(res$Var1[which(res$Term == cdr_info$Term[i])]))
cdr_info <- cdr_info[order(cdr_info$SN, decreasing=T),]

res$Term <- factor(res$Term, levels=cdr_info$Term, labels=paste0("(", cdr_info$SN, ")", cdr_info$Term))
res$Var1 <- factor(res$Var1, levels=sample_list)
resc <- res

pc <- ggplot(resc, aes(x=Term, y=Var1, size=Freq))+geom_point(color=col_list[19])+
	labs(title=NULL, x=NULL, y=NULL, size="Num. of\nBPlasma")+
	theme(plot.title=element_text(size=title_size, hjust=0.5, face="bold"), 
	axis.text.x=element_text(size=text_size*0.8, angle=270, hjust=0, vjust=0.5, colour="black"), 
	axis.text.y=element_text(size=text_size, colour="black"), panel.background=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	strip.background=element_rect(fill="white", colour="white"), 
	axis.ticks=element_blank(), axis.line=element_line(colour="black"))

bcrinfo_sub <- bcrinfo[which(bcrinfo$type == "BMemory" | bcrinfo$type == "BPlasma"),]
res_type <- data.frame()
for (term in cdr_info$Term) res_type <- rbind(res_type, data.frame(Term=term, 
	BM=length(unique(bcrinfo_sub$barcode[which(bcrinfo_sub$cdr3 == term & bcrinfo_sub$type == "BMemory")])), 
	BP=length(unique(bcrinfo_sub$barcode[which(bcrinfo_sub$cdr3 == term & bcrinfo_sub$type == "BPlasma")]))))
pa <- ggplot(res_type, aes(x=Term, y=BM))+geom_bar(stat="identity", width=0.8, position=position_dodge(0.8), fill=col_list[18])+
	labs(title=NULL, x=NULL, y=NULL)+ylim(0, max(res_type$BP, res_type$BM))+
	theme(axis.text.x=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	panel.background=element_blank(), axis.ticks.x=element_blank(), 
	axis.ticks.y=element_line(colour="black"), axis.line=element_line(colour="black"))
pb <- ggplot(res_type, aes(x=Term, y=BP))+geom_bar(stat="identity", width=0.8, position=position_dodge(0.8), fill=col_list[19])+
	labs(title=NULL, x=NULL, y=NULL)+ylim(0, max(res_type$BP, res_type$BM))+
	theme(axis.text.x=element_blank(), axis.text.y=element_text(size=text_size, colour="black"), 
	panel.background=element_blank(), axis.ticks.x=element_blank(), 
	axis.ticks.y=element_line(colour="black"), axis.line=element_line(colour="black"))

pd <- ggplot(data.frame(X=1:2, Y=1:2, Type=c("BMemory", "BPlasma")), aes(x=X, y=Y, fill=Type))+geom_bar(stat="identity")+
	scale_fill_manual(values=col_list[c(18,19)])+labs(title=NULL, fill="Type", x=NULL, y=NULL)+
	theme(panel.background=element_blank(), axis.text=element_blank(), 
	legend.title=element_text(size=title_size, colour="black"), 
	legend.text=element_text(size=text_size, colour="black"), 
	axis.ticks=element_blank(), axis.line=element_blank())

ggsave(plot=wrap_plots(list(pd, pa, pb, pc), ncol=1, heights=c(0,1,1,4))+plot_layout(guides="collect")+
	tag_thm+plot_annotation(tag_levels=list(c("I"))), 
	width=18, height=8, dpi=200, filename="cdr_test.png", limitsize=F)
ggsave(plot=wrap_plots(list(pd, pa, pb, pc), ncol=1, heights=c(0,1,1,4))+plot_layout(guides="collect")+
	tag_thm+plot_annotation(tag_levels=list(c("I"))), 
	width=18, height=8, dpi=200, filename="cdr_test.pdf", limitsize=F)

intersect(cdr_info$Term, c("CQQSYSPYTF", "CQQYYSTPLTF"))


ggsave(plot=wrap_plots(list(pd, pa, pb, pc), ncol=1, heights=c(0,1,1,4))+plot_layout(guides="collect"), 
	width=18, height=8, dpi=200, filename="cdr_test.png", limitsize=F)
