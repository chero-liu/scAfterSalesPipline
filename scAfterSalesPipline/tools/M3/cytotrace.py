source /home/liuxuan/miniconda3/bin/activate /gpfs/oe-software/conda_envs/scrna_envs/cytoTRACE
 
Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M3/cytotrace/CytoTrace_v1.2.R \
    -i /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023122645/20241031/data/Epithelial/Manualanno/data_ob_v3.rds   \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023122645/result/M3/CytoTrace/Paneth_2_C_1_Paneth_2_C_2_paneth_1_Stem_cell \
    -g new_celltype \
    -b FALSE \
    --subnew_celltype Paneth_2_C_1,Paneth_2_C_2,paneth_1,Stem_cell
