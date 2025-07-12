
#!/bin/bash
module purge && module load OESingleCell/2.0.0

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M4/infercnv//infercnv.R \
    --input asjhdkajhs.rds \
    --informat seurat \
    --celltype new_celltype \
    --refgroup ref \
    --gtf gtf \
    --clusting2use  all   \
    --output /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M4/inferCNV/Main


Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M4/infercnv//infercnv_vis.R \
    -i singlecell_object.clustering_resolution0.4.rds \
    -f  seurat \
    -l ./ \
    -g cnv_group,new_celltype,clusters,sampleid,group  \
    -p all \
    --reduct umap \
    -y sampleid,group \
    --hmm F

