
#!/bin/bash
module purge && module load OESingleCell/2.0.0

# Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/cellcycle/cellcyle_typing.R \
#     -m scran \
#     -i /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M1/Fibroblast_of_Main_20241127/Clustering/data_ob_v3.rds \
#     -f seurat \
#     -s mouse \
#     -b clusters,sampleid,group \
#     -t stack \
#     -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M3/CellCycle/Fibroblast_of_Main_20241127/ \
#     --reduct umap



Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/cellcycle/Cellcycle_3phases.r  \
    -i  /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M1/Fibroblast_of_Main_20241127/Clustering/data_ob_v3.rds  \
    -r umap \
    -s 0.5 \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M3/CellCycle/Seurat/Fibroblast_of_Main_20241127/ \
    --subgroup 

