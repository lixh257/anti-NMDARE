#an example 
cd /md01/zhangwx/project/human_PBMC/lib4/
/md01/yangjw28/software/cellranger-arc-2.0.1/bin/cellranger-arc count \
--id=lib4 \
--reference=/md01/zhangwx/ref/human/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries=/md01/zhangwx/rawdata/human_PBMC/lib4/libraries_lib4.csv \
--localcores=6 \
--localmem=80
