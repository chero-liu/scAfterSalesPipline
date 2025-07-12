module purge & module load OESingleCell/2.0.0

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/extra/homologene_transformed.R\
    -v singlecell_object.clustering_resolution0.4.rds \
    --assay RNA \
    -o ./ \
    -i 10116 \
    -t 10090


#!/bin/bash
module purge && module load OESingleCell/2.0.0

Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/homologene_transformed.R \
    -v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023112066/B3_Ileum/result/M1/Main/Manualanno/downsample_1000/data_ob_v3.rds \
    --assay RNA \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023112066/B3_Ileum/result/M1/Main/Manualanno/downsample_1000 \
    --blast T \
    --ingenome  "/gpfs/oe-database/cellranger-refdata/refdata-ARS-UI_Ramb_v2.0/fasta/genome.fa","/gpfs/oe-database/cellranger-refdata/refdata-ARS-UI_Ramb_v2.0/genes/genes.gtf.gz",/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023112066/B1/snakemake/DZOE2023112066-b1-Sheep/result/cellranger/D10JC/outs/filtered_feature_bc_matrix/features.tsv.gz \
    --outgenome /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/fasta/genome.fa,/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf,/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023050491/snakemake/result/cellranger/SG4/outs/filtered_feature_bc_matrix/features.tsv.gz


####Capra_hircus山羊
module purge && module load OESingleCell/2.0.0

Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/homologene_transformed.R \
    -v /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024091155/result/M1/Main/Manualanno/data_ob_v3.rds \
    --assay RNA \
    -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024091155/result/M1/Main/Manualanno \
    --blast T \
    --ingenome  "/gpfs/oe-database/cellranger-refdata/pre-mRNA/refdata-Capra_hircus-pre-mRNA/fasta/genome.fa","/gpfs/oe-database/cellranger-refdata/pre-mRNA/refdata-Capra_hircus-pre-mRNA/genes/genes.gtf",/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024091155/cellranger/s3_12/outs/filtered_feature_bc_matrix/features.tsv.gz \
    --outgenome /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/fasta/genome.fa,/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf,/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2023050491/snakemake/result/cellranger/SG4/outs/filtered_feature_bc_matrix/features.tsv.gz

