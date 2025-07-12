module purge & module load OESingleCell/3.0.d


for i in group:S:C
do

y=$(echo "$i" | awk -F: '{print $2 "_" $3}')
z=$(echo "$i" | awk -F: '{print $1}')
t=$(echo "$i" | awk -F: '{print $2 ":" $3}')
u=$(echo "$i" | awk -F: '{print $2 "," $3}')

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M5/CellChat_v1.6.1.R \
  -i /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024092169/2024-12-06/data/data_ob_v3.rds \
  -f rds \
  -s mouse  \
  -c new_celltype \
  -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024092169/result/M5/CellChat/Microglia_map_Main/$z/$y   \
  -g $z  \
  -d $t  \
  -q  $z \
  -u  $u

done
