module purge & module load OESingleCell/3.0.d


Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/CellChat_v1.6.1.R \
  -i /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024010167/result/M1/Main/Manualanno/data_ob_v3.rds \
  -f rds \
  -s mouse  \
  -c new_celltype \
  -o /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024010167/result/M5/cellchat/Main/CS-vs-RA   \
  -g group  \
  -d CS:RA  \
  -q group \
  -u CS,RA
