
#!/bin/bash
module purge && module load OESingleCell/2.0.0
Rscript   /home/luyao/10X_scRNAseq_v3/src/Diffexp/visualize_markers.R  \
-v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/result/M1/macrophage/Manualanno/data_ob_v3.rds \
-x /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/20240821/data/mic_cosscatter_infla.txt \
-o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/result/M2/Cosscatter/macrophage/mic_cosscatter_infla \
-m corscatter  \
-s 0.5 \
-g new_celltype \
-e RNA

# Rscript   /home/luyao/10X_scRNAseq_v3/src/Diffexp/visualize_markers.R  \
# -v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/result/M1/macrophage/Manualanno/data_ob_v3.rds \
# -x /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/20240821/data/m1m2.xls \
# -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024060483/result/M2/Cosscatter/macrophage/mic_cosscatter_m1m2 \
# -m corscatter  \
# -s 0.5 \
# -g new_celltype \
# -e RNA

