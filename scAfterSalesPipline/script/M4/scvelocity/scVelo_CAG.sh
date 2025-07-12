
module purge && module load velocity/1.0.0
velocyto run10x  /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023101502/20240712/data/exo_1 "/gpfs/oe-database/cellranger-refdata/refdata-gex-mm10-2020-A/genes/genes.gtf"

# /gpfs/oe-database/cellranger-refdata/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
module load git
module load scvelo
/gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M4/scvelocity/sctool \
    -i xx.rds \
    -f rds \
    --assay RNA \
    -j 10 \
    -o ./output  pyscvelo \
    --loom_dir yourloomdir \
    --groupby clusters,new_celltype \
    --reduction umap

