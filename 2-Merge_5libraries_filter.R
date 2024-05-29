#cellranger output file path:
Lib1:/data/R03/zhangwx/project/human_PBMC/lib1/outs
Lib2:/data/R03/zhangwx/project/human_PBMC/lib2/lib2/outs
Lib3:/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs
Lib4:/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs
Lib5:/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs


##in this part
1.filter low quality cells form every Library(filter doublet by scDblFinder)
2. merge
3.filter doublet and split samples by Mitosort(exclude Lib2)
4 Add sex of each samples



1.filter every Library
#scATAC:https://www.sc-best-practices.org/chromatin_accessibility/quality_control.html
#total_fragment_counts: Total number of fragments per cell representing cellular sequencing depth. /nCount_ATAC
#tss_enrichment: Transcription start site (TSS) enrichment score, which is the ratio of fragments centered at the TSS to fragments in TSS-flanking regions. 
#n_features_per_cell: The number of peaks with non-zero counts in each cell.
#nucleosome_signal: The nucleosome signal refers to the ratio of mono-nucleosomal to nucloesome-free fragments and can also be interpreted as a signal-to-noise ratio in each cell.
#Additional metrics that can be considered:
#reads_in_peaks_frac: The fraction of fragments in peak regions versus fragments outside of peaks. Similar to the TSS score, this is an indicator for the signal-to-noise ratio.
#blacklist_fraction: The ratio of fragments in genomic blacklist regions, which have been associated with artefactual signal 

##scRNA:https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
#The number of counts per barcode (count depth)
#The number of genes per barcode
#The fraction of counts from mitochondrial genes per barcode

##last cutoff in this study
#atac_fragments[50000,300]
#TSS.enrichment > 3
#nucleosome_signal < 2
#FRiP>0.2     
#nCount_RNA[5000 ,500] 
#gex_genes_count [3000, 200]
#percent.mt < 10 




1.1 Lib1
#data source location
#/data/R03/zhangwx/project/human_PBMC/lib1/outs/
1.1.1 Creat seurat project
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb) 
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) ##GRCh38
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(hdf5r)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1")
inputdata.lib1<- Read10X_h5(filename = "/data/R03/zhangwx/project/human_PBMC/lib1/outs/filtered_feature_bc_matrix.h5") 
rna_counts <- inputdata.lib1$`Gene Expression`
atac_counts <- inputdata.lib1$Peaks
metadata <- read.csv(file = "/data/R03/zhangwx/project/human_PBMC/lib1/outs/per_barcode_metrics.csv",   header = TRUE,   row.names = 1 )
# Create Seurat object 
pbmc_lib1 <- CreateSeuratObject(counts = rna_counts,meta.data = metadata) 
pbmc_lib1[["percent.mt"]] <- PercentageFeatureSet(pbmc_lib1, pattern = "^MT-")
head(pbmc_lib1@meta.data)

# Now add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) 
atac_counts <- atac_counts[as.vector(grange.use), ] 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotations) <- 'UCSC' 
genome(annotations) <- "GRCh38"

frag.file <- "/data/R03/zhangwx/project/human_PBMC/lib1/outs/atac_fragments.tsv.gz" 
chrom_assay <- CreateChromatinAssay( 
counts = atac_counts, 
sep = c(":", "-"), 
genome = 'GRCh38', 
fragments = frag.file, 
min.cells = 10, 
annotation = annotations
 ) 
pbmc_lib1[["ATAC"]] <- chrom_assay

# call peaks using MACS2
DefaultAssay(pbmc_lib1) <- "ATAC" 
peaks <- CallPeaks(pbmc_lib1, macs2.path = "/md01/liyh526/miniconda3/bin/macs2")
macs2_counts <- FeatureMatrix(
      fragments = Fragments(pbmc_lib1),
      features = peaks,
      cells = colnames(pbmc_lib1)
    )
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc_lib1[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = frag.file,
      annotation = annotations
    )

#pbmc_lib1[['peaks']]
DefaultAssay(pbmc_lib1) <- "peaks" 
granges(pbmc_lib1)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC' # add the gene information to the object
Annotation(pbmc_lib1) <- annotations

pbmc_lib1 <- NucleosomeSignal(object =pbmc_lib1) # compute TSS enrichment score per cell 
pbmc_lib1 <- TSSEnrichment(object =pbmc_lib1, fast = TRUE) # add blacklist ratio and fraction of reads in peaks 

#add FRiP
pbmc_lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
Lib1_total_fragments <- CountFragments("/data/R03/zhangwx/project/human_PBMC/lib1/outs/atac_fragments.tsv.gz")
rownames(Lib1_total_fragments) <- Lib1_total_fragments$CB
pbmc_lib1$fragments <- Lib1_total_fragments[colnames(pbmc_lib1), "frequency_count"]
pbmc_lib1 <- FRiP(
  object = pbmc_lib1,
  assay = 'peaks',
  total.fragments = 'fragments'
)
saveRDS(pbmc_lib1,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")

###filter doublet
###所以先将过低的cell filter,不然会影响doublet计算
library(scDblFinder)
pbmc_lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
pbmc_lib1_filter1<- subset( 
x = pbmc_lib1, 
atac_fragments > 300 &
nCount_RNA > 500 & 
gex_genes_count > 200 
)

###需要先进行降维分类等生成metadata
pbmc_lib1=pbmc_lib1_filter
DefaultAssay(pbmc_lib1) <- "RNA"
pbmc_lib1 <- SCTransform(pbmc_lib1, verbose = FALSE)
#integrate RNA using rpca
pbmc_lib1 <- FindVariableFeatures(pbmc_lib1)
pbmc_lib1 <- ScaleData(pbmc_lib1)
pbmc_lib1.sce.rna <- as.SingleCellExperiment(pbmc_lib1) 
pbmc_lib1 <- RunPCA(pbmc_lib1)
pbmc_lib1 <- RunUMAP(pbmc_lib1, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(pbmc_lib1) <- "ATAC"
pbmc_lib1 <- RunTFIDF(pbmc_lib1)
pbmc_lib1.sce.atac <- as.SingleCellExperiment(pbmc_lib1) 
pbmc_lib1 <- FindTopFeatures(pbmc_lib1, min.cutoff = "q25")
pbmc_lib1 <- RunSVD(pbmc_lib1)
pbmc_lib1 <- RunUMAP(pbmc_lib1, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_lib1 <- FindMultiModalNeighbors(pbmc_lib1, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc_lib1 <- RunUMAP(pbmc_lib1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_lib1 <- FindClusters(pbmc_lib1, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(pbmc_lib1$seurat_clusters)

###用RNA计算doublet
##预测的doublet比例大概是dbr <- weighted.mean(length(cell)/100000, length(cell))
pbmc_lib1.sce.rna <- scDblFinder(pbmc_lib1.sce.rna,cluster=pbmc_lib1$seurat_clusters)#, dbr=0.1
table(pbmc_lib1.sce.rna$scDblFinder.class)

pbmc_lib1.sce.atac <- scDblFinder(pbmc_lib1.sce.atac,cluster=pbmc_lib1$seurat_clusters)#, dbr=0.17364
table(pbmc_lib1.sce.atac$scDblFinder.class)

data=cbind(pbmc_lib1.sce.rna$scDblFinder.class,pbmc_lib1.sce.rna$scDblFinder.weighted,
           pbmc_lib1.sce.rna$scDblFinder.score,
           pbmc_lib1.sce.atac$scDblFinder.class,pbmc_lib1.sce.atac$scDblFinder.weighted,
           pbmc_lib1.sce.atac$scDblFinder.score)
colnames(data)=c("doublet","weight","score","atac_doublet","atac_weight","atac_score")
data=as.data.frame(data)
data[which(data[,1]==1),1]="singlet"
data[which(data[,1]==2),1]="doublet"
data[which(data[,4]==1),4]="singlet"
data[which(data[,4]==2),4]="doublet"

pbmc_lib1@meta.data$scDblFinder.rna=rep("Singlet",length(pbmc_lib1@meta.data$TSS.enrichment))
pbmc_lib1@meta.data$scDblFinder.rna[which(pbmc_lib1.sce.rna$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib1@meta.data$scDblFinder.atac=rep("Singlet",length(pbmc_lib1@meta.data$TSS.enrichment))
pbmc_lib1@meta.data$scDblFinder.atac[which(pbmc_lib1.sce.atac$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib1@meta.data$scDblFinder=rep("Singlet",length(pbmc_lib1@meta.data$TSS.enrichment))
pbmc_lib1@meta.data$scDblFinder[which(pbmc_lib1$scDblFinder.rna=="Doublet"|pbmc_lib1$scDblFinder.atac=="Doublet")]="Doublet"
doublet_barcode=rownames(pbmc_lib1@meta.data)[which(pbmc_lib1@meta.data$scDblFinder=="Doublet")]

###保持数据是原始状态
pbmc_lib1_raw=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
pbmc_lib1_raw@meta.data$scDblFinder.rna.atac=rep("Singlet",length(pbmc_lib1_raw@meta.data$TSS.enrichment))
pbmc_lib1_raw@meta.data$scDblFinder.rna.atac[which(rownames(pbmc_lib1_raw@meta.data)%in%doublet_barcode)]="Doublet"
saveRDS(pbmc_lib1_raw,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")

1.1.2 QC
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/")
pbmc_lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
pdf("pbmc_lib1_QC.pdf",width=15,height=6)
VlnPlot(   object = pbmc_lib1,   features = c('atac_fragments',"FRiP", 'TSS.enrichment', 
  'nucleosome_signal','percent.mt','nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'),
   pt.size = 0,   ncol = 9 )
dev.off()

cell_num=length(pbmc_lib1$TSS.enrichment)
mean_TSS=mean(pbmc_lib1@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib1@meta.data$TSS.enrichment)
fragment_max=max(pbmc_lib1@meta.data$atac_fragments)
fragment_min=min(pbmc_lib1@meta.data$atac_fragments)
fragment_median=median(pbmc_lib1@meta.data$atac_fragments)
Mean_raw_reads_GEX=mean(pbmc_lib1@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib1@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib1@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib1_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)
write.table(QC_table,"Lib1_QC_table.txt",quote=F,row.names=F,col.names=F,sep="\t")


pdf("Lib1_cutoff.pdf",width=10,height=4)
par(mfrow=c(1,4))
plot(density(pbmc_lib1@meta.data$nCount_RNA),xlim=c(0,5000))
plot(density(pbmc_lib1@meta.data$gex_genes_count),xlim=c(0,3000))
plot(density(pbmc_lib1@meta.data$nFeature_ATAC),xlim=c(0,15000))
plot(density(pbmc_lib1@meta.data$atac_fragments),xlim=c(0,50000))
dev.off()

1.1.3 filter 
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/")
pbmc_lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1.rds")
pbmc_lib1
pbmc_lib1_filter<- subset( 
x = pbmc_lib1, 
atac_fragments < 50000 &
atac_fragments > 300 &
TSS.enrichment > 3 &
nucleosome_signal < 2 &
FRiP>0.2 &      
nCount_RNA < 5000 & 
nCount_RNA > 500 & 
gex_genes_count < 3000 &
gex_genes_count > 200 &
percent.mt < 10 &
scDblFinder.rna.atac!="Doublet"
)
pbmc_lib1_filter

cell_num=length(pbmc_lib1_filter$TSS.enrichment)
fragment_max=max(pbmc_lib1_filter@meta.data$atac_fragments)
fragment_min=min(pbmc_lib1_filter@meta.data$atac_fragments)
fragment_median=median(pbmc_lib1_filter@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib1_filter@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib1_filter@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib1_filter@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib1_filter@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib1_filter@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib1_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)

write.table(QC_table,"Lib1_QC_table_filter.txt",quote=F,row.names=F,col.names=F,sep="\t")
saveRDS(pbmc_lib1_filter,"lib1_filter.rds")


1.2 Lib2
#Data: /md01/zhangwx/project/human_PBMC/lib2/lib2/outs/
1.2.1 
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) ##GRCh38
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(hdf5r)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2")
inputdata.lib2<- Read10X_h5(filename = "/md01/zhangwx/project/human_PBMC/lib2/lib2/outs/filtered_feature_bc_matrix.h5") 
rna_counts <- inputdata.lib2$`Gene Expression`
atac_counts <- inputdata.lib2$Peaks
metadata <- read.csv(file = "/md01/zhangwx/project/human_PBMC/lib2/lib2/outs/per_barcode_metrics.csv",   header = TRUE,   row.names = 1 )
# Create Seurat object 
pbmc_lib2 <- CreateSeuratObject(counts = rna_counts,meta.data = metadata) 
pbmc_lib2[["percent.mt"]] <- PercentageFeatureSet(pbmc_lib2, pattern = "^MT-")
head(pbmc_lib2@meta.data)

# Now add in the ATAC-seq data 
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) 
atac_counts <- atac_counts[as.vector(grange.use), ] 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotations) <- 'UCSC' 
genome(annotations) <- "GRCh38"

frag.file <- "/md01/zhangwx/project/human_PBMC/lib2/lib2/outs/atac_fragments.tsv.gz" 
chrom_assay <- CreateChromatinAssay( 
counts = atac_counts, 
sep = c(":", "-"), 
genome = 'GRCh38', 
fragments = frag.file, 
min.cells = 10, 
annotation = annotations
 ) 
pbmc_lib2[["ATAC"]] <- chrom_assay
 
# call peaks using MACS2
DefaultAssay(pbmc_lib2) <- "ATAC" 
peaks <- CallPeaks(pbmc_lib2, macs2.path = "/md01/liyh526/miniconda3/bin/macs2")
macs2_counts <- FeatureMatrix(
      fragments = Fragments(pbmc_lib2),
      features = peaks,
      cells = colnames(pbmc_lib2)
    )
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc_lib2[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = frag.file,
      annotation = annotations
    )

#pbmc_lib1[['peaks']]
DefaultAssay(pbmc_lib2) <- "peaks" 
granges(pbmc_lib2)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # change to UCSC style since the data was mapped to hg19
 seqlevelsStyle(annotations) <- 'UCSC' # add the gene information to the object
 Annotation(pbmc_lib2) <- annotations

pbmc_lib2 <- NucleosomeSignal(object =pbmc_lib2) 
pbmc_lib2 <- TSSEnrichment(object =pbmc_lib2, fast = TRUE) # # compute TSS enrichment score per cell 
#add blacklist ratio and fraction of reads in peaks 
#add FRiP
pbmc_lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")
Lib2_total_fragments <- CountFragments("/md01/zhangwx/project/human_PBMC/lib2/lib2/outs/atac_fragments.tsv.gz")
rownames(Lib2_total_fragments) <- Lib2_total_fragments$CB
pbmc_lib2$fragments <- Lib2_total_fragments[colnames(pbmc_lib2), "frequency_count"]
pbmc_lib2 <- FRiP(
  object = pbmc_lib2,
  assay = 'peaks',
  total.fragments = 'fragments'
)
saveRDS(pbmc_lib2,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")


######这个Libiary没有pooling ,只有一个样本，所以没办法使用MitSort的预测Doublet
##Warning message:
#In .checkSCE(sce) :
#  Some cells in `sce` have an extremely low read counts; note that these could trigger errors and might best be filtered out
###所以先将过低的cell filter
###filter doublet
###所以先将过低的cell filter,不然会影响doublet计算
library(scDblFinder)
pbmc_lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")
pbmc_lib2_filter1<- subset( 
x = pbmc_lib2, 
atac_fragments > 300 &
nCount_RNA > 500 & 
gex_genes_count > 200 
)

###需要先进行降维分类等生成metadata
pbmc_lib2=pbmc_lib2_filter1
DefaultAssay(pbmc_lib2) <- "RNA"
pbmc_lib2 <- SCTransform(pbmc_lib2, verbose = FALSE)
#integrate RNA using rpca
pbmc_lib2 <- FindVariableFeatures(pbmc_lib2)
pbmc_lib2 <- ScaleData(pbmc_lib2)
pbmc_lib2.sce.rna <- as.SingleCellExperiment(pbmc_lib2) 
pbmc_lib2 <- RunPCA(pbmc_lib2)
pbmc_lib2 <- RunUMAP(pbmc_lib2, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(pbmc_lib2) <- "ATAC"
pbmc_lib2 <- RunTFIDF(pbmc_lib2)
pbmc_lib2.sce.atac <- as.SingleCellExperiment(pbmc_lib2) 
pbmc_lib2 <- FindTopFeatures(pbmc_lib2, min.cutoff = "q25")
pbmc_lib2 <- RunSVD(pbmc_lib2)
pbmc_lib2 <- RunUMAP(pbmc_lib2, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_lib2 <- FindMultiModalNeighbors(pbmc_lib2, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc_lib2 <- RunUMAP(pbmc_lib2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_lib2 <- FindClusters(pbmc_lib2, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(pbmc_lib2$seurat_clusters)

###用RNA计算doublet
##预测的doublet比例大概是dbr <- weighted.mean(length(cell)/100000, length(cell))
pbmc_lib2.sce.rna <- scDblFinder(pbmc_lib2.sce.rna,cluster=pbmc_lib2$seurat_clusters)#, dbr=0.1
table(pbmc_lib2.sce.rna$scDblFinder.class)

pbmc_lib2.sce.atac <- scDblFinder(pbmc_lib2.sce.atac,cluster=pbmc_lib2$seurat_clusters)#, dbr=0.17364
table(pbmc_lib2.sce.atac$scDblFinder.class)

data=cbind(pbmc_lib2.sce.rna$scDblFinder.class,pbmc_lib2.sce.rna$scDblFinder.weighted,
           pbmc_lib2.sce.rna$scDblFinder.score,
           pbmc_lib2.sce.atac$scDblFinder.class,pbmc_lib2.sce.atac$scDblFinder.weighted,
           pbmc_lib2.sce.atac$scDblFinder.score)
colnames(data)=c("doublet","weight","score","atac_doublet","atac_weight","atac_score")
data=as.data.frame(data)
data[which(data[,1]==1),1]="singlet"
data[which(data[,1]==2),1]="doublet"
data[which(data[,4]==1),4]="singlet"
data[which(data[,4]==2),4]="doublet"
write.table(data,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/scDblFinder_result.txt",sep="\t",quote=F)

pbmc_lib2@meta.data$scDblFinder.rna=rep("Singlet",length(pbmc_lib2@meta.data$TSS.enrichment))
pbmc_lib2@meta.data$scDblFinder.rna[which(pbmc_lib2.sce.rna$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib2@meta.data$scDblFinder.atac=rep("Singlet",length(pbmc_lib2@meta.data$TSS.enrichment))
pbmc_lib2@meta.data$scDblFinder.atac[which(pbmc_lib2.sce.atac$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib2@meta.data$scDblFinder=rep("Singlet",length(pbmc_lib2@meta.data$TSS.enrichment))
pbmc_lib2@meta.data$scDblFinder[which(pbmc_lib2$scDblFinder.rna=="Doublet"|pbmc_lib2$scDblFinder.atac=="Doublet")]="Doublet"
doublet_barcode=rownames(pbmc_lib2@meta.data)[which(pbmc_lib2@meta.data$scDblFinder=="Doublet")]
length(doublet_barcode)
###保持数据是原始状态
pbmc_lib2_raw=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")
pbmc_lib2_raw@meta.data$scDblFinder.rna.atac=rep("Singlet",length(pbmc_lib2_raw@meta.data$TSS.enrichment))
pbmc_lib2_raw@meta.data$scDblFinder.rna.atac[which(rownames(pbmc_lib2_raw@meta.data)%in%doublet_barcode)]="Doublet"
table(pbmc_lib2_raw$scDblFinder.rna.atac)
saveRDS(pbmc_lib2_raw,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")


1.2.2 QC
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/")
pbmc_lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")
pdf("pbmc_lib2_QC.pdf",width=15,height=6)
VlnPlot(   object = pbmc_lib2,   features = c('atac_fragments',"FRiP", 'TSS.enrichment', 
  'nucleosome_signal','percent.mt','nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'),
   pt.size = 0,   ncol = 9 )
dev.off()

cell_num=length(pbmc_lib2$TSS.enrichment)
mean_TSS=mean(pbmc_lib2@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib2@meta.data$TSS.enrichment)
fragment_max=max(pbmc_lib2@meta.data$atac_fragments)
fragment_min=min(pbmc_lib2@meta.data$atac_fragments)
fragment_median=median(pbmc_lib2@meta.data$atac_fragments)
Mean_raw_reads_GEX=mean(pbmc_lib2@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib2@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib2@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib2_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)
write.table(QC_table,"Lib2_QC_table.txt",quote=F,row.names=F,col.names=F,sep="\t")

#
pdf("Lib2_cutoff.pdf",width=10,height=4)
par(mfrow=c(1,4))
plot(density(pbmc_lib2@meta.data$nCount_RNA),xlim=c(0,5000))
plot(density(pbmc_lib2@meta.data$gex_genes_count),xlim=c(0,3000))
plot(density(pbmc_lib2@meta.data$nFeature_ATAC),xlim=c(0,15000))
plot(density(pbmc_lib2@meta.data$atac_fragments),xlim=c(0,50000))
dev.off()

1.2.3 filter
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/")
pbmc_lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2.rds")

##filter doublet
#pbmc_lib2=subset(pbmc_lib2,scDblFinder.rna.atac!="Doublet")
pbmc_lib2_filter<- subset( 
x = pbmc_lib2, 
atac_fragments < 50000 &
atac_fragments > 300 &
TSS.enrichment > 3 &
nucleosome_signal < 2 &
FRiP>0.2 &      
nCount_RNA < 5000 & 
nCount_RNA > 500 & 
gex_genes_count < 3000 &
gex_genes_count > 200 &
percent.mt < 10 &
scDblFinder.rna.atac!="Doublet"
)
pbmc_lib2_filter

cell_num=length(pbmc_lib2_filter$TSS.enrichment)
fragment_max=max(pbmc_lib2_filter@meta.data$atac_fragments)
fragment_min=min(pbmc_lib2_filter@meta.data$atac_fragments)
fragment_median=median(pbmc_lib2_filter@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib2_filter@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib2_filter@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib2_filter@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib2_filter@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib2_filter@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib2_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)

write.table(QC_table,"Lib2_QC_table_filter.txt",quote=F,row.names=F,col.names=F,sep="\t")
saveRDS(pbmc_lib2_filter,"lib2_filter.rds")



1.3 Lib3
1.3 Lib3_Add
#Source:/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs
1.3.1
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb) 
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) ##GRCh38
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(hdf5r)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/")
inputdata.lib3<- Read10X_h5(filename = "/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/filtered_feature_bc_matrix.h5") 
rna_counts <- inputdata.lib3$`Gene Expression`
atac_counts <- inputdata.lib3$Peaks
metadata <- read.csv(file = "/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/per_barcode_metrics.csv",   header = TRUE,   row.names = 1 )
# Create Seurat object 
pbmc_lib3 <- CreateSeuratObject(counts = rna_counts,meta.data = metadata) 
pbmc_lib3[["percent.mt"]] <- PercentageFeatureSet(pbmc_lib3, pattern = "^MT-")
head(pbmc_lib3@meta.data)

# Now add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) 
atac_counts <- atac_counts[as.vector(grange.use), ] 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotations) <- 'UCSC' 
genome(annotations) <- "GRCh38"

frag.file <- "/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/atac_fragments.tsv.gz" 
chrom_assay <- CreateChromatinAssay( 
counts = atac_counts, 
sep = c(":", "-"), 
genome = 'GRCh38', 
fragments = frag.file, 
min.cells = 10, 
annotation = annotations
 ) 
pbmc_lib3[["ATAC"]] <- chrom_assay

# call peaks using MACS2
DefaultAssay(pbmc_lib3) <- "ATAC" 
peaks <- CallPeaks(pbmc_lib3, macs2.path = "/md01/liyh526/miniconda3/bin/macs2")
macs2_counts <- FeatureMatrix(
      fragments = Fragments(pbmc_lib3),
      features = peaks,
      cells = colnames(pbmc_lib3)
    )
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc_lib3[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = frag.file,
      annotation = annotations
    )

#pbmc_lib1[['peaks']]
DefaultAssay(pbmc_lib3) <- "peaks" 
granges(pbmc_lib3)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC' # add the gene information to the object
Annotation(pbmc_lib3) <- annotations

pbmc_lib3 <- NucleosomeSignal(object =pbmc_lib3) # compute TSS enrichment score per cell 
pbmc_lib3 <- TSSEnrichment(object =pbmc_lib3, fast = TRUE) # add blacklist ratio and fraction of reads in peaks 
#add FRiP
Lib3_total_fragments <- CountFragments("/data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/atac_fragments.tsv.gz")
rownames(Lib3_total_fragments) <- Lib3_total_fragments$CB
pbmc_lib3$fragments <- Lib3_total_fragments[colnames(pbmc_lib3), "frequency_count"]
pbmc_lib3 <- FRiP(
  object = pbmc_lib3,
  assay = 'peaks',
  total.fragments = 'fragments'
)
saveRDS(pbmc_lib3,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")


###filter doublet
###所以先将过低的cell filter,不然会影响doublet计算
library(scDblFinder)
pbmc_lib3=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")
pbmc_lib3_filter1<- subset( 
x = pbmc_lib3, 
atac_fragments > 300 &
nCount_RNA > 500 & 
gex_genes_count > 200 
)

###需要先进行降维分类等生成metadata
pbmc_lib3=pbmc_lib3_filter1
DefaultAssay(pbmc_lib3) <- "RNA"
pbmc_lib3 <- SCTransform(pbmc_lib3, verbose = FALSE)
#integrate RNA using rpca
pbmc_lib3 <- FindVariableFeatures(pbmc_lib3)
pbmc_lib3 <- ScaleData(pbmc_lib3)
pbmc_lib3.sce.rna <- as.SingleCellExperiment(pbmc_lib3) 
pbmc_lib3 <- RunPCA(pbmc_lib3)
pbmc_lib3 <- RunUMAP(pbmc_lib3, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(pbmc_lib3) <- "ATAC"
pbmc_lib3 <- RunTFIDF(pbmc_lib3)
pbmc_lib3.sce.atac <- as.SingleCellExperiment(pbmc_lib3) 
pbmc_lib3 <- FindTopFeatures(pbmc_lib3, min.cutoff = "q25")
pbmc_lib3 <- RunSVD(pbmc_lib3)
pbmc_lib3 <- RunUMAP(pbmc_lib3, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_lib3 <- FindMultiModalNeighbors(pbmc_lib3, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc_lib3 <- RunUMAP(pbmc_lib3, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_lib3 <- FindClusters(pbmc_lib3, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(pbmc_lib3$seurat_clusters)

###用RNA计算doublet
##预测的doublet比例大概是dbr <- weighted.mean(length(cell)/100000, length(cell))
pbmc_lib3.sce.rna <- scDblFinder(pbmc_lib3.sce.rna,cluster=pbmc_lib3$seurat_clusters)#, dbr=0.1
table(pbmc_lib3.sce.rna$scDblFinder.class)

pbmc_lib3.sce.atac <- scDblFinder(pbmc_lib3.sce.atac,cluster=pbmc_lib3$seurat_clusters)#, dbr=0.17364
table(pbmc_lib3.sce.atac$scDblFinder.class)

data=cbind(pbmc_lib3.sce.rna$scDblFinder.class,pbmc_lib3.sce.rna$scDblFinder.weighted,
           pbmc_lib3.sce.rna$scDblFinder.score,
           pbmc_lib3.sce.atac$scDblFinder.class,pbmc_lib3.sce.atac$scDblFinder.weighted,
           pbmc_lib3.sce.atac$scDblFinder.score)
colnames(data)=c("doublet","weight","score","atac_doublet","atac_weight","atac_score")
data=as.data.frame(data)
data[which(data[,1]==1),1]="singlet"
data[which(data[,1]==2),1]="doublet"
data[which(data[,4]==1),4]="singlet"
data[which(data[,4]==2),4]="doublet"
write.table(data,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/scDblFinder_result.txt",sep="\t",quote=F)

pbmc_lib3@meta.data$scDblFinder.rna=rep("Singlet",length(pbmc_lib3@meta.data$TSS.enrichment))
pbmc_lib3@meta.data$scDblFinder.rna[which(pbmc_lib3.sce.rna$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib3@meta.data$scDblFinder.atac=rep("Singlet",length(pbmc_lib3@meta.data$TSS.enrichment))
pbmc_lib3@meta.data$scDblFinder.atac[which(pbmc_lib3.sce.atac$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib3@meta.data$scDblFinder=rep("Singlet",length(pbmc_lib3@meta.data$TSS.enrichment))
pbmc_lib3@meta.data$scDblFinder[which(pbmc_lib3$scDblFinder.rna=="Doublet"|pbmc_lib3$scDblFinder.atac=="Doublet")]="Doublet"
doublet_barcode=rownames(pbmc_lib3@meta.data)[which(pbmc_lib3@meta.data$scDblFinder=="Doublet")]
length(doublet_barcode)
###保持数据是原始状态
pbmc_lib3_raw=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")
pbmc_lib3_raw@meta.data$scDblFinder.rna.atac=rep("Singlet",length(pbmc_lib3_raw@meta.data$TSS.enrichment))
pbmc_lib3_raw@meta.data$scDblFinder.rna.atac[which(rownames(pbmc_lib3_raw@meta.data)%in%doublet_barcode)]="Doublet"
table(pbmc_lib3_raw$scDblFinder.rna.atac)
saveRDS(pbmc_lib3_raw,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")

1.3.2 QC
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/")
pbmc_lib3=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3.rds")
pdf("pbmc_lib3_QC.pdf",width=15,height=6)
VlnPlot(   object = pbmc_lib3,   features = c('atac_fragments',"FRiP", 'TSS.enrichment', 
  'nucleosome_signal','percent.mt','nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'),
   pt.size = 0,   ncol = 9 )
dev.off()


cell_num=length(pbmc_lib3$TSS.enrichment)
mean_TSS=mean(pbmc_lib3@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib3@meta.data$TSS.enrichment)
fragment_max=max(pbmc_lib3@meta.data$atac_fragments)
fragment_min=min(pbmc_lib3@meta.data$atac_fragments)
fragment_median=median(pbmc_lib3@meta.data$atac_fragments)
Mean_raw_reads_GEX=mean(pbmc_lib3@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib3@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib3@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib3_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)
write.table(QC_table,"Lib3_QC_table.txt",quote=F,row.names=F,col.names=F,sep="\t")


##determination the cutoff
pdf("Lib3_cutoff.pdf",width=10,height=4)
par(mfrow=c(1,4))
plot(density(pbmc_lib3@meta.data$nCount_RNA),xlim=c(0,5000))
plot(density(pbmc_lib3@meta.data$gex_genes_count),xlim=c(0,3000))
plot(density(pbmc_lib3@meta.data$nFeature_ATAC),xlim=c(0,15000))
plot(density(pbmc_lib3@meta.data$atac_fragments),xlim=c(0,50000))
dev.off()

pbmc_lib3
pbmc_lib3_filter<- subset( 
x = pbmc_lib3, 
atac_fragments < 50000 &
atac_fragments > 300 &
TSS.enrichment > 3 &
nucleosome_signal < 2 &
FRiP>0.2 &      
nCount_RNA < 5000 & 
nCount_RNA > 500 & 
gex_genes_count < 3000 &
gex_genes_count > 200 &
percent.mt < 10 &
scDblFinder.rna.atac!="Doublet"
)
pbmc_lib3_filter

cell_num=length(pbmc_lib3_filter$TSS.enrichment)
fragment_max=max(pbmc_lib3_filter@meta.data$atac_fragments)
fragment_min=min(pbmc_lib3_filter@meta.data$atac_fragments)
fragment_median=median(pbmc_lib3_filter@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib3_filter@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib3_filter@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib3_filter@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib3_filter@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib3_filter@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib3_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)

write.table(QC_table,"Lib3_QC_table_filter.txt",quote=F,row.names=F,col.names=F,sep="\t")
saveRDS(pbmc_lib3_filter,"lib3_filter.rds")




1.4 Lib4_Add
#source:/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs
1.4.1
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb) 
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) ##GRCh38
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(hdf5r)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/")
inputdata.lib4<- Read10X_h5(filename = "/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/filtered_feature_bc_matrix.h5") 
rna_counts <- inputdata.lib4$`Gene Expression`
atac_counts <- inputdata.lib4$Peaks
metadata <- read.csv(file = "/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/per_barcode_metrics.csv",   header = TRUE,   row.names = 1 )
# Create Seurat object 
pbmc_lib4 <- CreateSeuratObject(counts = rna_counts,meta.data = metadata) 
pbmc_lib4[["percent.mt"]] <- PercentageFeatureSet(pbmc_lib4, pattern = "^MT-")
head(pbmc_lib4@meta.data)

# Now add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) 
atac_counts <- atac_counts[as.vector(grange.use), ] 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotations) <- 'UCSC' 
genome(annotations) <- "GRCh38"

frag.file <- "/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/atac_fragments.tsv.gz" 
chrom_assay <- CreateChromatinAssay( 
counts = atac_counts, 
sep = c(":", "-"), 
genome = 'GRCh38', 
fragments = frag.file, 
min.cells = 10, 
annotation = annotations
 ) 
pbmc_lib4[["ATAC"]] <- chrom_assay

# call peaks using MACS2
DefaultAssay(pbmc_lib4) <- "ATAC" 
peaks <- CallPeaks(pbmc_lib4, macs2.path = "/md01/liyh526/miniconda3/bin/macs2")
macs2_counts <- FeatureMatrix(
      fragments = Fragments(pbmc_lib4),
      features = peaks,
      cells = colnames(pbmc_lib4)
    )
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc_lib4[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = frag.file,
      annotation = annotations
    )

#pbmc_lib1[['peaks']]
DefaultAssay(pbmc_lib4) <- "peaks" 
granges(pbmc_lib4)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC' # add the gene information to the object
Annotation(pbmc_lib4) <- annotations

pbmc_lib4 <- NucleosomeSignal(object =pbmc_lib4) # compute TSS enrichment score per cell 
pbmc_lib4 <- TSSEnrichment(object =pbmc_lib4, fast = TRUE) # add blacklist ratio and fraction of reads in peaks 

##add FRiP
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/")
pbmc_lib4=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")
Lib4_total_fragments <- CountFragments("/data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/atac_fragments.tsv.gz")
rownames(Lib4_total_fragments) <- Lib4_total_fragments$CB
pbmc_lib4$fragments <- Lib4_total_fragments[colnames(pbmc_lib4), "frequency_count"]
pbmc_lib4 <- FRiP(
  object = pbmc_lib4,
  assay = 'peaks',
  total.fragments = 'fragments'
)
saveRDS(pbmc_lib4,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")

###filter doublet
###所以先将过低的cell filter,不然会影响doublet计算
library(scDblFinder)
pbmc_lib4=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")
pbmc_lib4_filter1<- subset( 
x = pbmc_lib4, 
atac_fragments > 300 &
nCount_RNA > 500 & 
gex_genes_count > 200 
)

###需要先进行降维分类等生成metadata
pbmc_lib4=pbmc_lib4_filter1
DefaultAssay(pbmc_lib4) <- "RNA"
pbmc_lib4 <- SCTransform(pbmc_lib4, verbose = FALSE)
#integrate RNA using rpca
pbmc_lib4 <- FindVariableFeatures(pbmc_lib4)
pbmc_lib4 <- ScaleData(pbmc_lib4)
pbmc_lib4.sce.rna <- as.SingleCellExperiment(pbmc_lib4) 
pbmc_lib4 <- RunPCA(pbmc_lib4)
pbmc_lib4 <- RunUMAP(pbmc_lib4, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(pbmc_lib4) <- "ATAC"
pbmc_lib4 <- RunTFIDF(pbmc_lib4)
pbmc_lib4.sce.atac <- as.SingleCellExperiment(pbmc_lib4) 
pbmc_lib4 <- FindTopFeatures(pbmc_lib4, min.cutoff = "q25")
pbmc_lib4 <- RunSVD(pbmc_lib4)
pbmc_lib4 <- RunUMAP(pbmc_lib4, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_lib4 <- FindMultiModalNeighbors(pbmc_lib4, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc_lib4 <- RunUMAP(pbmc_lib4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_lib4 <- FindClusters(pbmc_lib4, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(pbmc_lib4$seurat_clusters)

###用RNA计算doublet
##预测的doublet比例大概是dbr <- weighted.mean(length(cell)/100000, length(cell))
pbmc_lib4.sce.rna <- scDblFinder(pbmc_lib4.sce.rna,cluster=pbmc_lib4$seurat_clusters)#, dbr=0.1
table(pbmc_lib4.sce.rna$scDblFinder.class)

pbmc_lib4.sce.atac <- scDblFinder(pbmc_lib4.sce.atac,cluster=pbmc_lib4$seurat_clusters)#, dbr=0.17364
table(pbmc_lib4.sce.atac$scDblFinder.class)

data=cbind(pbmc_lib4.sce.rna$scDblFinder.class,pbmc_lib4.sce.rna$scDblFinder.weighted,
           pbmc_lib4.sce.rna$scDblFinder.score,
           pbmc_lib4.sce.atac$scDblFinder.class,pbmc_lib4.sce.atac$scDblFinder.weighted,
           pbmc_lib4.sce.atac$scDblFinder.score)
colnames(data)=c("doublet","weight","score","atac_doublet","atac_weight","atac_score")
data=as.data.frame(data)
data[which(data[,1]==1),1]="singlet"
data[which(data[,1]==2),1]="doublet"
data[which(data[,4]==1),4]="singlet"
data[which(data[,4]==2),4]="doublet"

pbmc_lib4@meta.data$scDblFinder.rna=rep("Singlet",length(pbmc_lib4@meta.data$TSS.enrichment))
pbmc_lib4@meta.data$scDblFinder.rna[which(pbmc_lib4.sce.rna$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib4@meta.data$scDblFinder.atac=rep("Singlet",length(pbmc_lib4@meta.data$TSS.enrichment))
pbmc_lib4@meta.data$scDblFinder.atac[which(pbmc_lib4.sce.atac$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib4@meta.data$scDblFinder=rep("Singlet",length(pbmc_lib4@meta.data$TSS.enrichment))
pbmc_lib4@meta.data$scDblFinder[which(pbmc_lib4$scDblFinder.rna=="Doublet"|pbmc_lib4$scDblFinder.atac=="Doublet")]="Doublet"
doublet_barcode=rownames(pbmc_lib4@meta.data)[which(pbmc_lib4@meta.data$scDblFinder=="Doublet")]

###保持数据是原始状态
pbmc_lib4_raw=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")
pbmc_lib4_raw@meta.data$scDblFinder.rna.atac=rep("Singlet",length(pbmc_lib4_raw@meta.data$TSS.enrichment))
pbmc_lib4_raw@meta.data$scDblFinder.rna.atac[which(rownames(pbmc_lib4_raw@meta.data)%in%doublet_barcode)]="Doublet"
saveRDS(pbmc_lib4_raw,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")


1.4.2 QC
pbmc_lib4=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4.rds")
pdf("pbmc_lib4_QC.pdf",width=15,height=6)
VlnPlot(   object = pbmc_lib4,   features = c('atac_fragments',"FRiP", 'TSS.enrichment', 
  'nucleosome_signal','percent.mt','nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'),
   pt.size = 0,   ncol = 9 )
dev.off()

cell_num=length(pbmc_lib4$TSS.enrichment)
mean_TSS=mean(pbmc_lib4@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib4@meta.data$TSS.enrichment)
fragment_max=max(pbmc_lib4@meta.data$atac_fragments)
fragment_min=min(pbmc_lib4@meta.data$atac_fragments)
fragment_median=median(pbmc_lib4@meta.data$atac_fragments)
Mean_raw_reads_GEX=mean(pbmc_lib4@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib4@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib4@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib4_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)
write.table(QC_table,"Lib4_QC_table.txt",quote=F,row.names=F,col.names=F,sep="\t")

1.4.3 filter
##determination the cutoff
pdf("Lib4_cutoff.pdf",width=10,height=4)
par(mfrow=c(1,4))
plot(density(pbmc_lib4@meta.data$nCount_RNA),xlim=c(0,5000))
plot(density(pbmc_lib4@meta.data$gex_genes_count),xlim=c(0,3000))
plot(density(pbmc_lib4@meta.data$nFeature_ATAC),xlim=c(0,15000))
plot(density(pbmc_lib4@meta.data$atac_fragments),xlim=c(0,50000))
dev.off()

pbmc_lib4
pbmc_lib4_filter<- subset( 
x = pbmc_lib4, 
atac_fragments < 50000 &
atac_fragments > 300 &
TSS.enrichment > 3 &
nucleosome_signal < 2 &
FRiP>0.2 &      
nCount_RNA < 5000 & 
nCount_RNA > 500 & 
gex_genes_count < 3000 &
gex_genes_count > 200 &
percent.mt < 10 &
scDblFinder.rna.atac!="Doublet"
)
pbmc_lib4_filter

cell_num=length(pbmc_lib4_filter$TSS.enrichment)
fragment_max=max(pbmc_lib4_filter@meta.data$atac_fragments)
fragment_min=min(pbmc_lib4_filter@meta.data$atac_fragments)
fragment_median=median(pbmc_lib4_filter@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib4_filter@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib4_filter@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib4_filter@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib4_filter@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib4_filter@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib4_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)

write.table(QC_table,"Lib4_QC_table_filter.txt",quote=F,row.names=F,col.names=F,sep="\t")
saveRDS(pbmc_lib4_filter,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4_filter.rds")


1.5 Lib5
/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/
1.5.1
library(Seurat)
library(Signac)
library(dplyr) 
library(GenomeInfoDb) 
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) ##GRCh38
library(ggplot2)
library(Hmisc)
library(biovizBase)
library(lmtest)
#library(rsvg)
library(hdf5r)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5")
 inputdata.lib5<- Read10X_h5(filename = "/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/filtered_feature_bc_matrix.h5") 
rna_counts <- inputdata.lib5$`Gene Expression`
 atac_counts <- inputdata.lib5$Peaks
metadata <- read.csv(file = "/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/per_barcode_metrics.csv",   header = TRUE,   row.names = 1 )
# Create Seurat object 
pbmc_lib5 <- CreateSeuratObject(counts = rna_counts,meta.data = metadata) 
pbmc_lib5[["percent.mt"]] <- PercentageFeatureSet(pbmc_lib5, pattern = "^MT-")
head(pbmc_lib5@meta.data)

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes 
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) 
atac_counts <- atac_counts[as.vector(grange.use), ] 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevelsStyle(annotations) <- 'UCSC' 
genome(annotations) <- "GRCh38"

frag.file <- "/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/atac_fragments.tsv.gz" 
chrom_assay <- CreateChromatinAssay( 
counts = atac_counts, 
sep = c(":", "-"), 
genome = 'GRCh38', 
fragments = frag.file, 
min.cells = 10, 
annotation = annotations
 ) 
pbmc_lib5[["ATAC"]] <- chrom_assay

# call peaks using MACS2
DefaultAssay(pbmc_lib5) <- "ATAC" 
peaks <- CallPeaks(pbmc_lib5, macs2.path = "/md01/liyh526/miniconda3/bin/macs2")
macs2_counts <- FeatureMatrix(
      fragments = Fragments(pbmc_lib5),
      features = peaks,
      cells = colnames(pbmc_lib5)
    )
# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc_lib5[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = frag.file,
      annotation = annotations
    )

#pbmc_lib1[['peaks']]
DefaultAssay(pbmc_lib5) <- "peaks" 
granges(pbmc_lib5)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # change to UCSC style since the data was mapped to hg19
 seqlevelsStyle(annotations) <- 'UCSC' # add the gene information to the object
 Annotation(pbmc_lib5) <- annotations

pbmc_lib5 <- NucleosomeSignal(object =pbmc_lib5) # compute TSS enrichment score per cell 
pbmc_lib5 <- TSSEnrichment(object =pbmc_lib5, fast = TRUE) # add blacklist ratio and fraction of reads in peaks 

pbmc_lib5=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")
Lib5_total_fragments <- CountFragments("/data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/atac_fragments.tsv.gz")
rownames(Lib5_total_fragments) <- Lib5_total_fragments$CB
pbmc_lib5$fragments <- Lib5_total_fragments[colnames(pbmc_lib5), "frequency_count"]
pbmc_lib5 <- FRiP(
  object = pbmc_lib5,
  assay = 'peaks',
  total.fragments = 'fragments'
)
saveRDS(pbmc_lib5,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")

###filter doublet
###所以先将过低的cell filter,不然会影响doublet计算
library(scDblFinder)
pbmc_lib5=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")
pbmc_lib5_filter1<- subset( 
x = pbmc_lib5, 
atac_fragments > 300 &
nCount_RNA > 500 & 
gex_genes_count > 200 
)

###需要先进行降维分类等生成metadata
pbmc_lib5=pbmc_lib5_filter1
DefaultAssay(pbmc_lib5) <- "RNA"
pbmc_lib5 <- SCTransform(pbmc_lib5, verbose = FALSE)
#integrate RNA using rpca
pbmc_lib5 <- FindVariableFeatures(pbmc_lib5)
pbmc_lib5 <- ScaleData(pbmc_lib5)
pbmc_lib5.sce.rna <- as.SingleCellExperiment(pbmc_lib5) 
pbmc_lib5 <- RunPCA(pbmc_lib5)
pbmc_lib5 <- RunUMAP(pbmc_lib5, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#run LSI on new seurat object with integrated RNA assay
DefaultAssay(pbmc_lib5) <- "ATAC"
pbmc_lib5 <- RunTFIDF(pbmc_lib5)
pbmc_lib5.sce.atac <- as.SingleCellExperiment(pbmc_lib5) 
pbmc_lib5 <- FindTopFeatures(pbmc_lib5, min.cutoff = "q25")
pbmc_lib5 <- RunSVD(pbmc_lib5)
pbmc_lib5 <- RunUMAP(pbmc_lib5, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc_lib5 <- FindMultiModalNeighbors(pbmc_lib5, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
pbmc_lib5 <- RunUMAP(pbmc_lib5, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_lib5 <- FindClusters(pbmc_lib5, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.5)
table(pbmc_lib5$seurat_clusters)

###用RNA计算doublet
##预测的doublet比例大概是dbr <- weighted.mean(length(cell)/100000, length(cell))
pbmc_lib5.sce.rna <- scDblFinder(pbmc_lib5.sce.rna,cluster=pbmc_lib5$seurat_clusters)#, dbr=0.1
table(pbmc_lib5.sce.rna$scDblFinder.class)

pbmc_lib5.sce.atac <- scDblFinder(pbmc_lib5.sce.atac,cluster=pbmc_lib5$seurat_clusters)#, dbr=0.17364
table(pbmc_lib5.sce.atac$scDblFinder.class)

data=cbind(pbmc_lib5.sce.rna$scDblFinder.class,pbmc_lib5.sce.rna$scDblFinder.weighted,
           pbmc_lib5.sce.rna$scDblFinder.score,
           pbmc_lib5.sce.atac$scDblFinder.class,pbmc_lib5.sce.atac$scDblFinder.weighted,
           pbmc_lib5.sce.atac$scDblFinder.score)
colnames(data)=c("doublet","weight","score","atac_doublet","atac_weight","atac_score")
data=as.data.frame(data)
data[which(data[,1]==1),1]="singlet"
data[which(data[,1]==2),1]="doublet"
data[which(data[,4]==1),4]="singlet"
data[which(data[,4]==2),4]="doublet"

pbmc_lib5@meta.data$scDblFinder.rna=rep("Singlet",length(pbmc_lib5@meta.data$TSS.enrichment))
pbmc_lib5@meta.data$scDblFinder.rna[which(pbmc_lib5.sce.rna$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib5@meta.data$scDblFinder.atac=rep("Singlet",length(pbmc_lib5@meta.data$TSS.enrichment))
pbmc_lib5@meta.data$scDblFinder.atac[which(pbmc_lib5.sce.atac$scDblFinder.class=="doublet")]="Doublet"

pbmc_lib5@meta.data$scDblFinder=rep("Singlet",length(pbmc_lib5@meta.data$TSS.enrichment))
pbmc_lib5@meta.data$scDblFinder[which(pbmc_lib5$scDblFinder.rna=="Doublet"|pbmc_lib5$scDblFinder.atac=="Doublet")]="Doublet"
doublet_barcode=rownames(pbmc_lib5@meta.data)[which(pbmc_lib5@meta.data$scDblFinder=="Doublet")]

###保持数据是原始状态
pbmc_lib5_raw=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")
pbmc_lib5_raw@meta.data$scDblFinder.rna.atac=rep("Singlet",length(pbmc_lib5_raw@meta.data$TSS.enrichment))
pbmc_lib5_raw@meta.data$scDblFinder.rna.atac[which(rownames(pbmc_lib5_raw@meta.data)%in%doublet_barcode)]="Doublet"
table(pbmc_lib5_raw@meta.data$scDblFinder.rna.atac)
saveRDS(pbmc_lib5_raw,"/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")


1.5.2 QC
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/")
pbmc_lib5=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5.rds")
pdf("pbmc_lib5_QC.pdf",width=15,height=6)
VlnPlot(   object = pbmc_lib5,   features = c('atac_fragments',"FRiP", 'TSS.enrichment', 
  'nucleosome_signal','percent.mt','nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'),
   pt.size = 0,   ncol = 9 )
dev.off()


cell_num=length(pbmc_lib5$TSS.enrichment)
fragment_max=max(pbmc_lib5@meta.data$atac_fragments)
fragment_min=min(pbmc_lib5@meta.data$atac_fragments)
fragment_median=median(pbmc_lib5@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib5@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib5@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib5@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib5@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib5@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib5_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)
write.table(QC_table,"Lib5_QC_table.txt",quote=F,row.names=F,col.names=F,sep="\t")

1.5.3 filter
##determination the cutoff
pdf("Lib5_cutoff.pdf",width=10,height=4)
par(mfrow=c(1,4))
plot(density(pbmc_lib5@meta.data$nCount_RNA),xlim=c(0,5000))
plot(density(pbmc_lib5@meta.data$gex_genes_count),xlim=c(0,3000))
plot(density(pbmc_lib5@meta.data$nFeature_ATAC),xlim=c(0,15000))
plot(density(pbmc_lib5@meta.data$atac_fragments),xlim=c(0,50000))
dev.off()

pbmc_lib5
pbmc_lib5_filter<- subset( 
x = pbmc_lib5, 
atac_fragments < 50000 &
atac_fragments > 300 &
TSS.enrichment > 3 &
nucleosome_signal < 2 &
FRiP>0.2 &      
nCount_RNA < 5000 & 
nCount_RNA > 500 & 
gex_genes_count < 3000 &
gex_genes_count > 200 &
percent.mt < 10 &
scDblFinder.rna.atac!="Doublet"
)
pbmc_lib5_filter

cell_num=length(pbmc_lib5_filter$TSS.enrichment)
fragment_max=max(pbmc_lib5_filter@meta.data$atac_fragments)
fragment_min=min(pbmc_lib5_filter@meta.data$atac_fragments)
fragment_median=median(pbmc_lib5_filter@meta.data$atac_fragments)
mean_TSS=mean(pbmc_lib5_filter@meta.data$TSS.enrichment)
median_TSS=median(pbmc_lib5_filter@meta.data$TSS.enrichment)
Mean_raw_reads_GEX=mean(pbmc_lib5_filter@meta.data$gex_raw_reads)
Mean_raw_reads_ATAC=mean(pbmc_lib5_filter@meta.data$atac_raw_reads)
median_gex_genes_count=median(pbmc_lib5_filter@meta.data$gex_genes_count)
QC_table=matrix(0,ncol=2,nrow=8)
QC_table[1,]=c("Lib5_cell_number",cell_num)
QC_table[2,]=c("mean_TSS",mean_TSS)
QC_table[3,]=c("median_TSS",median_TSS)
QC_table[4,]=c("fragment_scope",paste0(fragment_min,"-",fragment_max))
QC_table[5,]=c("fragment_median",fragment_median)
QC_table[6,]=c("Mean_raw_reads_GEX",Mean_raw_reads_GEX)
QC_table[7,]=c("Mean_raw_reads_ATAC",Mean_raw_reads_ATAC)
QC_table[8,]=c("median_gex_genes_count",median_gex_genes_count)

write.table(QC_table,"Lib5_QC_table_filter.txt",quote=F,row.names=F,col.names=F,sep="\t")
saveRDS(pbmc_lib5_filter,"lib5_filter.rds")





2.merge
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
library(dittoSeq)
library(dplyr)
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
Lib1=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib1/lib1_filter.rds")
Lib2=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib2/lib2_filter.rds")
Lib3=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib3/Lib3_Add/lib3_filter.rds")
Lib4=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib4/Lib4_Add/lib4_filter.rds")
Lib5=readRDS("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Lib5/lib5_filter.rds")

Lib1$dataset <- 'Lib1' 
Lib2$dataset <- 'Lib2' 
Lib3$dataset <- 'Lib3' 
Lib4$dataset <- 'Lib4' 
Lib5$dataset <- 'Lib5'
Merge_after <- merge(   x = Lib1,   y = list(Lib2,Lib3,Lib4,Lib5),   add.cell.ids = c("Lib1","Lib2", "Lib3", "Lib4","Lib5") )
saveRDS(Merge_after,"Merge_after_filter1.rds")



3.filter doublet and split samples by Mitosort
3.1 MitoSort
Lib1:
python ./MitoSort_pipeline.py demultiplex -o \
/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib1 \
-k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 --depth_cutoff 1 --method 'direct'


Lib3:
python ./MitoSort_pipeline.py demultiplex -o \
/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib3 \
-k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 --depth_cutoff 1 --method 'direct'

#########
Lib4:
python ./MitoSort_pipeline.py demultiplex -o \
/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib4 \
-k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 --depth_cutoff 1 --method 'direct'
#########

#########
Lib5:
python ./MitoSort_pipeline.py demultiplex -o \
/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/MitoSort_test/lib5 \
-k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 --depth_cutoff 1 --method 'direct'


3.2 split samples and delet doublets
#copy MitoSort result to my path"
#source: /data/R04/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/ 
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
combined=readRDS("Merge_after_filter1.rds")
pre_Lib1=read.table("/data/R04/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/lib1/MitoSort/Demultiplex_output/lib1_depth_1_result_pvalue.txt",sep="\t",head=T)
pre_Lib1_patient1=pre_Lib1[which(pre_Lib1[,2]=="Sample0"),1]
pre_Lib1_patient2=pre_Lib1[which(pre_Lib1[,2]=="Sample1"),1]

pre_Lib3=read.table("/data/R04/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/lib3/MitoSort/Demultiplex_output/lib3_depth_1_result_pvalue.txt",sep="\t",head=T)
pre_Lib3_patient=pre_Lib3[which(pre_Lib3[,2]=="Sample1"),1]#female,吴小丹
pre_Lib3_normal=pre_Lib3[which(pre_Lib3[,2]=="Sample0"),1]#male,Normal

pre_Lib4=read.table("/data/R04/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/lib4/MitoSort/Demultiplex_output/lib4_depth_1_result_pvalue.txt",sep="\t",head=T)
pre_Lib4_patient=pre_Lib4[which(pre_Lib4[,2]=="Sample1"),1]#female,冯芷琳
pre_Lib4_normal=pre_Lib4[which(pre_Lib4[,2]=="Sample0"),1]#male,Normal

pre_Lib5=read.table("/data/R04/chenbzh5/bottleneck2/simulation/bam_data/simulate_depth/MitoSort_test/lib5/MitoSort/Demultiplex_output/lib5_depth_1_result_pvalue.txt",sep="\t",head=T)
pre_Lib5_normal1=pre_Lib5[which(pre_Lib5[,2]=="Sample2"),1]
pre_Lib5_normal2=pre_Lib5[which(pre_Lib5[,2]=="Sample1"),1]
pre_Lib5_normal3=pre_Lib5[which(pre_Lib5[,2]=="Sample3"),1]
pre_Lib5_normal4=pre_Lib5[which(pre_Lib5[,2]=="Sample0"),1]

combined$Label=0
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib1"),which(substring(rownames(combined@meta.data),6,24) %in% pre_Lib1_patient1))]="Patient1"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib1"),which(substring(rownames(combined@meta.data),6,24) %in% pre_Lib1_patient2))]="Patient2"
combined$Label[which(substring(rownames(combined@meta.data),1,4)=="Lib2")]="Patient3"##ZCH
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib3"),which(substring(rownames(combined@meta.data),6,24) %in% pre_Lib3_patient))]="Patient4"##WXD
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib3"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib3_normal))]="Normal1"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib4"),which(substring(rownames(combined@meta.data),6,24)  %in%  pre_Lib4_patient))]="Patient5"##FZL
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib4"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib4_normal))]="Normal2"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib5"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib5_normal1))]="Normal3"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib5"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib5_normal2))]="Normal4"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib5"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib5_normal3))]="Normal5"
combined$Label[intersect(which(substring(rownames(combined@meta.data),1,4)=="Lib5"),which(substring(rownames(combined@meta.data),6,24)   %in% pre_Lib5_normal4))]="Normal6"
#> table(combined$Label)
#       0  Normal1  Normal2  Normal3  Normal4  Normal5  Normal6 Patient1 
#    4275     5327     5881     3221     4255     2870     3854    12491 
#Patient2 Patient3 Patient4 Patient5 
#     335     9775     6750     2987
###有4275个细胞不好分类或者是包含两个样本被删掉了


combined=subset(x=combined,Label != 0)
combined$AE=0
combined$AE[which(combined$Label %in%c("Patient1","Patient2","Patient3","Patient4","Patient5"))]="Patient"
combined$AE[which(combined$Label %in% c("Normal1","Normal2","Normal3","Normal4","Normal5","Normal6"))]="Normal"
table(combined$AE)


4 Add sex of samples
4.1 get chrX and chrY genes
cp -r /md01/zhangwx/ref/human/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz /md01/lixh/Autoimmune/Merge_4samples/
cd /md01/lixh/Autoimmune/Merge_4samples/
gunzip genes.gtf.gz 
grep "chrX" genes.gtf |perl -ne '/gene_name "(.*)"; transcript_type/;print"$1\n"'|sort|uniq > chrX.csv 
grep "chrY" genes.gtf |perl -ne '/gene_name "(.*)"; transcript_type/;print"$1\n"'|sort|uniq > chrY.csv

cp -r /md01/lixh/Autoimmune/Merge_4samples/genes.gtf /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/chrXY/
cp -r /md01/lixh/Autoimmune/Merge_4samples/chrX.csv /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/chrXY/
cp -r /md01/lixh/Autoimmune/Merge_4samples/chrY.csv /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/chrXY/

4.2 Add chrX and chrY score to combined
cY <- read.csv("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/chrXY/chrY.csv") 
#head(cY$gene) 
DefaultAssay(combined) <- "RNA" 
table(cY$gene %in% rownames(combined)) 

combined[['chrY']] <- PercentageFeatureSet(combined,features = cY$gene[cY$gene %in% rownames(combined)]) 
#table(combined1@meta.data$chrY>0.1) 

cX <- read.csv("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/chrXY/chrX.csv") 
head(cX$gene) 
#DefaultAssay(combined1) <- "RNA" 
table(cX$gene %in% rownames(combined)) 

combined[['chrX']] <- PercentageFeatureSet(combined,features = cX$gene[cX$gene %in% rownames(combined)]) 
table(combined@meta.data$chrX>0.1)
colnames(combined@meta.data)
chrXY=combined@meta.data[,c(46,48,49,50)]##dataset,Label,chrY,chrX
chrXY_Lib1=chrXY[which(chrXY[,1]=="Lib1"),]
chrXY_Lib1_P1=chrXY[which(chrXY[,2]=="Patient1"),]
chrXY_Lib1_P2=chrXY[which(chrXY[,2]=="Patient2"),]
table(chrXY_Lib1_P1[,3]>0)#Female
#FALSE  TRUE 
#12267   224
table(chrXY_Lib1_P2[,3]>0)#Male
#FALSE  TRUE 
#   34   301

chrXY_Lib2=chrXY[which(chrXY[,1]=="Lib2"),]
chrXY_Lib2_P=chrXY[which(chrXY[,2]=="Patient3"),]
table(chrXY_Lib2_P[,3]>0)#Male,也验证了可以用这个方法判断每个样本的性别
#FALSE  TRUE 
#  711  9064 


chrXY_Lib3=chrXY[which(chrXY[,1]=="Lib3"),]
chrXY_Lib3_P4=chrXY[which(chrXY[,2]=="Patient4"),]
chrXY_Lib3_N1=chrXY[which(chrXY[,2]=="Normal1"),]
table(chrXY_Lib3_P4[,3]>0)#Female，已知病人和样本的性别，所以借此判断哪个是病人哪个是健康对照
table(chrXY_Lib3_N1[,3]>0)#Male
#FALSE  TRUE 
# 6399   351 
#FALSE  TRUE 
# 1332  3995

chrXY_Lib4=chrXY[which(chrXY[,1]=="Lib4"),]
chrXY_Lib4_P5=chrXY[which(chrXY[,2]=="Patient5"),]
chrXY_Lib4_N2=chrXY[which(chrXY[,2]=="Normal2"),]
table(chrXY_Lib4_P5[,3]>0)#Female
#FALSE  TRUE 
# 2768   219 
table(chrXY_Lib4_N2[,3]>0)#Male
#FALSE  TRUE 
#  617  5264

chrXY_Lib5=chrXY[which(chrXY[,1]=="Lib5"),]
chrXY_Lib5_N3=chrXY[which(chrXY[,2]=="Normal3"),]
chrXY_Lib5_N4=chrXY[which(chrXY[,2]=="Normal4"),]
chrXY_Lib5_N5=chrXY[which(chrXY[,2]=="Normal5"),]
chrXY_Lib5_N6=chrXY[which(chrXY[,2]=="Normal6"),]
table(chrXY_Lib5_N3[,3]>0)#Female
table(chrXY_Lib5_N4[,3]>0)#Female
table(chrXY_Lib5_N5[,3]>0)#Male
table(chrXY_Lib5_N6[,3]>0)#Female


#combined=combined1
combined$Sex=0
combined$Sex[which(combined$Label=="Patient1")]="Female"#HYJ
combined$Sex[which(combined$Label=="Patient2")]="Male"#FSP
combined$Sex[which(combined$Label=="Patient3")]="Male"#ZCH
combined$Sex[which(combined$Label=="Patient4")]="Female"#WXD
combined$Sex[which(combined$Label=="Patient5")]="Female"#FZL
combined$Sex[which(combined$Label=="Normal1")]="Male"
combined$Sex[which(combined$Label=="Normal2")]="Male"
combined$Sex[which(combined$Label=="Normal3")]="Female"
combined$Sex[which(combined$Label=="Normal4")]="Female"
combined$Sex[which(combined$Label=="Normal5")]="Male"
combined$Sex[which(combined$Label=="Normal6")]="Female"
saveRDS(combined,"Merge_after_filter_includP2.rds")


combined=subset(x=combined,Label!="Patient2")##太少
saveRDS(combined,"Merge_after_filter1.rds")

###后期经过讨论决定要将Patient5（FZL，经过治疗）的样本去掉
##SO:
setwd("/data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/")
combined=readRDS("Merge_after_filter1.rds")
combined1=subset(combined,Label!="Patient5")
saveRDS(combined1,"Merge_after_filter1.rds")



