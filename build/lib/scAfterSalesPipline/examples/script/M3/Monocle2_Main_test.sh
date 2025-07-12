
#!/bin/bash
module purge && module load OESingleCell/2.0.0

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/monocle2/monocle.R \
    --INEXPRESS dksahdk.rds  \
    --INFORMAT seurat \
    --assay RNA \
    --design clusters \
    --CORES 8 \
    --min.cell 0.01 \
    --resolution 0.4 \
    --colorby clusters,sampleid,group,State \
    --pointsize 1 \
    --output /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M3/Monocle2/Main \
    --downsample 30000 \
    --batch  F


Rscript  /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/monocle2/visualize_pseudotime.R  \
-i /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2023122838/result/M3/monocle2/B_cells_vargene_bygroup/pseudotime_results.rds \
-m all \
-g ordering \
-c new_celltype \
--root  1  \
-b 2 \
--show_branch T \
-j 5 \
-a /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls \
-o /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2023122838/result/M3/monocle2/B_cells_vargene_bygroup/new_celltype


