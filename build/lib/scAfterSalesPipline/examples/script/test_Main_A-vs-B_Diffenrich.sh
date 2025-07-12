
#!/bin/bash
module purge && module load OESingleCell/3.0.d

file_path=/gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M1/Main/Diffexp/A-vs-B/done.check

while true; do
    if [ -f "$file_path" ]; then
        break
    else
        sleep 60
    fi
done

perl  /gpfs/oe-scrna/ziqingzhen/script/enrichment/enrich_go_kegg.pl  \
-infile  /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M1/Main/Diffexp/A-vs-B/group_*-diff-pval-0.05-FC-1.5.xls \
-go_bg /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_go.backgroud.xls \
-category  /home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls  \
-kegg_bg  /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_kegg.backgroud.xls  \
-outdir  /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/examples/result/M1/Main/Diffenrich  \
-shelldir enrichment_sh \
-thread 4 \
-queue big
