# for i in `find /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024020609/submission/DZOE2024020609-王宝祯-肺单细胞-外来数据分析_20240401 -type f -iname '*.pdf'`
# do
# echo $i
# python /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/extra/pdf_overlap_checker.py $i
# done
date
python /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/extra/pdf_overlap_checker.py /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024020609/submission/DZOE2024020609-王宝祯-肺单细胞-外来数据分析_20240401/T_NK_cell/Reference_celltype/humanref_hpca_single_celltyping_heatmap.pdf   --ignore_shape_overlap  --ignore_line           --ignore_curve          --ignore_rect           --ignore_quad           
date
