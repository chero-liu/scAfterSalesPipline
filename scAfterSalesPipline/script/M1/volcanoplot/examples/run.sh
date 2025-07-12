module purge && module load OESingleCell/3.0.d

Rscript ./volcano.r \
    -i examples/group_BAV-vs-AV-all_diffexp_genes.xls \
    -P BAV-vs-AV \
    -q 0.05 \
    -f 1.5 \
    -o examples/out_up_down \
    --symbol_topn 20


Rscript ./volcano.r \
    -i examples/group_BAV-vs-AV-all_diffexp_genes.xls \
    -P BAV-vs-AV \
    --symbol_gene  examples/BAV-vs-AV_show.xls \
    -q 0.05 \
    -f 1.5 \
    -o examples/out_label

