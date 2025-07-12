module load OESingleCell/3.0.d

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/extra/loupe.R \
    -i ./data_ob_v3.rds \
    -f rds \
    -o ./loupe \
    -p Mobi  # BD(BD项目)，Mobi(墨卓项目)，10× 不用填写-p
