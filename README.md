# anti-NMDARE(Anti-N-methyl-D-aspartate receptor encephalitis)
In this project, we processed the scRNA-seq,scATAC-seq and BCR-seq dataset of PBMCs form four anti-NMDARE patients and six healthy controls.

Parper :

Raw Data :

Code process:


1. Processing raw data by cellranger-arc-2.0.1
   
(1-cellranger-example.sh)


2.QC

(2-Merge_5libraries_filter.R)

##in this part:
Filter low quality cells form every Library(QC,filter doublet by scDblFinder);  
Merge 5 libraries data;  
Filter doublet and split samples by Mitosort(exclude Lib2);  
Add sex of each samples;  

3. Reduction and Annotation cell type

(3-Reduction_Annotation.R)

4-15. Figures

