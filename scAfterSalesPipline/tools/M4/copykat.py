Rscript  /home/luyao/10X_scRNAseq_v3/src/CNV/copykat.R  \
    -i  singlecell_object.clustering_resolution0.4.rds  \
    -l celltype  \
    --reduct umap \
    -g celltype \
    -o ./  \
    -j 10 -r T_cells
    