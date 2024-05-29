source activate Pyth3_R4
cd /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/test/

#lib1
#step1
python ./MitoSort_pipeline.py mt-realign -b /data/R03/zhangwx/project/human_PBMC/lib1/outs/atac_possorted_bam.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --gatk_path /md01/zhangwx/project/separeteLib1/MitoSort/GenomeAnalysisTK_3.5-0.jar -o /md01/zhangwx/project/separeteLib1
#step2
python ./MitoSort_pipeline.py generate-snp-matrix -b /md01/zhangwx/project/separeteLib1/MitoSort/BAM/possorted_chrM_realign.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa -c /data/R03/zhangwx/project/human_PBMC/lib1/outs/per_barcode_metrics.csv -m /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/data/hg38_chrM.bed --varscan_path /md01/zhangwx/project/separeteLib1/MitoSort/VarScan.v2.3.7.jar --cell_tag CB -o /md01/zhangwx/project/separeteLib1/
#step3
#routine pipeline for demultiplex sample
#python ./MitoSort_pipeline.py demultiplex -o /md01/zhangwx/project/separeteLib1/ -k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 

#because the minimal cell number of lib1 we should add the parmeter "--method 'direct'"
python ./MitoSort_pipeline.py demultiplex -o /md01/zhangwx/project/separeteLib1/ -k 2 --p1_cutoff 0.7 --p2_cutoff 0.3 --method 'direct'

#Lib3_add
#step1
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py mt-realign \
-b /data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/atac_possorted_bam.bam \
 -f /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
 --gatk_path /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/GenomeAnalysisTK_3.5-0.jar \
 -o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib3_add/

#step2
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib3_add/MitoSort/BAM/possorted_chrM.bam \
-f /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-c /data/R03/zhangwx/project/human_PBMC/lib3/joint_hPBMC_lib3_h38_add/outs/per_barcode_metrics.csv \
-m /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/data/hg38_chrM.bed \
--varscan_path /md01/zhangwx/project/separeteLib1/MitoSort/VarScan.v2.3.7.jar \
--cell_tag CB -o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib3_add/

#step3
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py demultiplex \
-o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib3_add -k 2 --p1_cutoff 0.8 --p2_cutoff 0.2

#Lib4_add
#step1
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py mt-realign \
-b /data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/atac_possorted_bam.bam \
 -f /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
 --gatk_path /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/GenomeAnalysisTK_3.5-0.jar \
 -o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib4_add/

#step2
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
-b /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib4_add/MitoSort/BAM/possorted_chrM.bam \
-f /data/R04/lixh/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
-c /data/R03/zhangwx/project/human_PBMC/lib4/lib4_add/lib4/outs/per_barcode_metrics.csv \
-m /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/data/hg38_chrM.bed \
--varscan_path /md01/zhangwx/project/separeteLib1/MitoSort/VarScan.v2.3.7.jar \
--cell_tag CB -o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib4_add/

#step3
python /data/R04/lixh/MitoSort/MitoSort_pipeline.py demultiplex \
-o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib4_add -k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 --method 'direct'

python /data/R04/lixh/MitoSort/MitoSort_pipeline.py demultiplex \
-o /data/R04/lixh/Autoimmune/1-5Lib_5AE_6Nor/Sup_Fig1_QC/Lib4_add -k 2 --p1_cutoff 1.0 --p2_cutoff 0


#lib5
#step1
python ./MitoSort_pipeline.py mt-realign 
-b /data/R03/zhangwx/project/human_PBMC/lib5/outs/atac_possorted_bam.bam
 -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa 
 --gatk_path /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/GenomeAnalysisTK_3.5-0.jar 
 o /md01/Chenzh275/Project/Zhangwx/lib5/

#step2
python ./MitoSort_pipeline.py generate-snp-matrix 
-b /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/BAM/possorted_chrM_realign.bam 
-f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa 
-c /data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/per_barcode_metrics.csv 
-m /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/data/hg38_chrM.bed 
--varscan_path /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/VarScan.v2.3.7.jar 
--cell_tag CB -o /md01/Chenzh275/Project/Zhangwx/lib5/

#step3
python ./MitoSort_pipeline.py demultiplex 
-o /md01/Chenzh275/Project/Zhangwx/lib5 -k 4 --p1_cutoff 0.8 --p2_cutoff 0.2


