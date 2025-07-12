#!/bin/bash
module purge && module load OESingleCell/2.0.0


Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M1/markerplot/sc_heatmap.r  \
    -i /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M1/B_Plasma_20240708/Clustering/data_ob_v3.rds  \
    -g /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/20240918/data/Bexpmarkerplot.txt  \
    -c clusters \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023113282/result/M1/B_Plasma_20240708/sc_heatmap/ \
    -r F \
    -l F \
    -a RNA \
    --topn 100 \
    --subcluster "2,3,4,6,8"



Rscript /home/luyao/10X_scRNAseq_v3/src/Diffexp/visualize_markers.R  \
    -v singlecell_object.clustering_resolution0.4.rds \
    -x genelist.txt \
    -o ./  \
    -m flip_vlnplot  \
    -s 0 \
    -e RNA \
    -g new_celltype \
    --reduct tsne
