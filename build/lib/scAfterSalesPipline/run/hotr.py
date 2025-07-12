module purge & module load OESingleCell/2.0.0

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/extra/homologene_transformed.R\
    -v singlecell_object.clustering_resolution0.4.rds \
    --assay RNA \
    -o ./ \
    -i 10116 \
    -t 10090
