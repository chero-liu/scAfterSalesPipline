
#!/bin/bash
module purge && module load OESingleCell/2.0.0

Rscript   /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M2/addmodule/visualize_markers.R  \
    --RDS akshdjkash.rds \
    --extraGene hskjahd.xls \
    --output /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M2/Addmodule/Main/hskjahd \
    --vismethod geneset  \
    --pointsize 0.5 \
    --assay RNA \
    --reduct umap  \
    --groupby clusters    

