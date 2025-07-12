module load OESingleCell/3.0.d

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/slingshot/slingshot.R \
    --RDS /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M1/Fibroblast_of_Main_20241127/Clustering/data_ob_v3.rds \
    --informat rds \
    --seurat F \
    --output /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024101066/result/M3/SlingShot/Fibroblast_of_Main_20241127/ \
    --groupby clusters \
    --start 1 \
    --pointsize 1 \
    --reduct UMAP
