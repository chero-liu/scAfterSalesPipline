import os

from scAfterSalesPipline.__init__ import ROOT_PATH, SOFTWARE_NAME

CONFIG = {
    "M1": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM1.yaml",
    "M2": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM2.yaml",
    "M3": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM3.yaml",
    "M4": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM4.yaml",
    "M5": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM5.yaml",
    "M6": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/templates/config/cfgM6.yaml",
}
SCRIPT = {
    "inferCNV": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M4/infercnv",
    "v4tov3": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/extra",
    "Monocle2": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M3/monocle2",
    "GSVA": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M2/gsva",
    "Addmodule": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M2/addmodule",
    "Clustering": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/clustering",
    "Manualanno": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/manualanno",
    "Markerplot": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/markerplot",
    "Diffexp": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/diffexp",
    "Enrich": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/enrich",
    "PPI": f"{os.path.dirname(ROOT_PATH)}/{SOFTWARE_NAME}/script/M1/ppi",
}
ENVIRONMENT = {
    "Clustering": "module purge && module load OESingleCell/3.0.d",
    "Manualanno": "module purge && module load OESingleCell/3.0.d",
    "Markerplot": "module purge && module load OESingleCell/3.0.d",
    "Diffexp": "module purge && module load OESingleCell/3.0.d",
    "Enrich": "module purge && module load OESingleCell/3.0.d",
    "PPI": "module purge && module load OESingleCell/3.0.d",
    "v4tov3": "module purge && module load OESingleCell/3.0.d",
    "GSVA": "module purge && module load OESingleCell/2.0.0",
    "Addmodule": "module purge && module load OESingleCell/2.0.0",
    "inferCNV": "module purge && module load OESingleCell/2.0.0",
    "Monocle2": "module purge && module load OESingleCell/2.0.0",
}

NCORES = {
    "Clustering": 5,
    "Manualanno": 5,
    "Markerplot": 2,
    "Diffexp": 3,
    "Enrich": 3,
    "PPI": 2,
    "v4tov3": 5,
    "GSVA": 20,
    "Addmodule": 3,
    "inferCNV": 20,
    "Monocle2": 20,
}

GENOME = {
    "go_bg": {
        "human": "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_go.backgroud.xls",
        "mouse": "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_go.backgroud.xls",
        "rat": "/gpfs/oe-database/cellranger-refdata/refdata-mRatBN7/annotation/gene_go.backgroud.xls",
        "gallus_gallus": "/data/database/cellranger-refdata/refdata-Gallus_gallus/GRCg7b/annotation/gene_go.backgroud.xls",
        "dog": "/data/database/cellranger-refdata/refdata-Canis_lupus_familiaris_v3.1/annotation/gene_go.backgroud.xls",
        "Caenorhabditis_elegans": "/gpfs/oe-scrna/liuchenglong/geneAnno/Caenorhabditis_elegans/annotation/gene_go.backgroud.xls",
        "oar": "/data/database/cellranger-refdata/refdata-ARS-UI_Ramb_v2.0/annotation/gene_go.backgroud.xls",
    },
    "kegg_bg": {
        "human": "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_kegg.backgroud.xls",
        "mouse": "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_kegg.backgroud.xls",
        "rat": "/gpfs/oe-database/cellranger-refdata/refdata-mRatBN7/annotation/gene_kegg.backgroud.xls",
        "gallus_gallus": "/data/database/cellranger-refdata/refdata-Gallus_gallus/GRCg7b/annotation/gene_kegg.backgroud.xls",
        "dog": "/data/database/cellranger-refdata/refdata-Canis_lupus_familiaris_v3.1/annotation/gene_kegg.backgroud.xls",
        "Caenorhabditis_elegans": "/gpfs/oe-scrna/liuchenglong/geneAnno/Caenorhabditis_elegans/annotation/gene_kegg.backgroud.xls",
        "oar": "/data/database/cellranger-refdata/refdata-ARS-UI_Ramb_v2.0/annotation/gene_kegg.backgroud.xls",
    },
    "geneanno": {
        "human": "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls",
        "mouse": "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_annotation.xls",
        "rat": "/gpfs/oe-database/cellranger-refdata/refdata-mRatBN7/annotation/gene_annotation.xls",
        "oar": "/gpfs/oe-database/cellranger-refdata/refdata-Oar_v4.0/annotation/gene_annotation.xls",
        "capra_hircus": "/data/database/cellranger-refdata/pre-mRNA/refdata-Capra_hircus-pre-mRNA/annotation/gene_annotation.xls",
        "mmul": "/data/database/cellranger-refdata/refdata-Mmul_10/annotation/gene_annotation.xls",
        "macaca_fascicularis": "/data/database/cellranger-refdata/refdata-Macaca_fascicularis_GCF_000364345.1_5.0/annotation/gene_annotation.xls",
        "sus_scrofa": "/data/database/cellranger-refdata/refdata-Sus_scrofa/annotation/gene_annotation.xls",
        "gallus_gallus": "/data/database/cellranger-refdata/refdata-Gallus_gallus/GRCg7b/annotation/gene_annotation.xls",
        "Caenorhabditis_elegans": "/gpfs/oe-scrna/liuchenglong/geneAnno/Caenorhabditis_elegans/annotation/gene_annotation.xls",
        "dog": "/data/database/cellranger-refdata/refdata-Canis_lupus_familiaris_v3.1/annotation/gene_annotation.xls",
    },
    "gtf": {
        "human": "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        "mouse": "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A/genes/genes.gtf",
        "rat": "/gpfs/oe-database/cellranger-refdata/refdata-mRatBN7/genes/genes.gtf",
        "oar": "/gpfs/oe-database/cellranger-refdata/refdata-Oar_v4.0/genes/genes.gtf",
        "capra_hircus": "/data/database/cellranger-refdata/pre-mRNA/refdata-Capra_hircus-pre-mRNA/genes/genes.gtf",
        "mmul": "/data/database/cellranger-refdata/refdata-Mmul_10/genes/genes.gtf",
        "macaca_fascicularis": "/data/database/cellranger-refdata/refdata-Macaca_fascicularis_GCF_000364345.1_5.0/genes/genes.gtf",
        "sus_scrofa": "/data/database/cellranger-refdata/refdata-Sus_scrofa/genes/genes.gtf",
        "gallus_gallus": "/data/database/cellranger-refdata/refdata-Gallus_gallus/GRCg7b/genes/genes.gtf",
        "Caenorhabditis_elegans": "/gpfs/oe-scrna/liuchenglong/geneAnno/Caenorhabditis_elegans/genes/genes.gtf",
        "dog": "/data/database/cellranger-refdata/refdata-Canis_lupus_familiaris_v3.1/genes/genes.gtf",
    },
    "cellanno": {
        "human": "/gpfs/oe-database/celltype_refdata/logNorm_rds/hpca.rds",
        "mouse": "/data/database/celltype_refdata/logNorm_rds/immgen.rds",
        "rat": "/data/database/celltype_refdata/logNorm_rds/immgen.rds",
        "oar": None,
        "capra_hircus": None,
        "mmul": "/gpfs/oe-database/celltype_refdata/logNorm_rds/hpca.rds",
        "macaca_fascicularis": "/gpfs/oe-database/celltype_refdata/logNorm_rds/hpca.rds",
        "sus_scrofa": "/gpfs/oe-database/celltype_refdata/logNorm_rds/hpca.rds",
        "gallus_gallus": "/gpfs/oe-database/celltype_refdata/logNorm_rds/hpca.rds",
        "dog": None,
        "Caenorhabditis_elegans": None,
    },
    "ppispecies": {
        "human": "9606",
        "mouse": "10090",
        "rat": "10116",
        "gallus_gallus": "9031",
        "dog": "9615",
        "Caenorhabditis_elegans": None,
        "oar": "9940",
    },
}


GMT = {
    "human": {
        "go": "/data/database/GSEA_gmt/human/v2023/c5.go.bp.v2023.1.Hs.symbols.gmt",
        "kegg": "/data/database/GSEA_gmt/human/v2023/c2.cp.kegg.v2023.1.Hs.symbols.gmt",
        "Hallmakr": "/public/scRNA_works/works/lipeng/script/hallmark/h.all.v2023.1.Hs.symbols.gmt",
    },
    "mouse": {
        "go": "/data/database/GSEA_gmt/mouse/v2023/m5.go.bp.v2023.1.Mm.symbols.gmt",
        "kegg": "/data/database/GSEA_gmt/mouse/gene_kegg.gmt",
        "Hallmakr": "/public/scRNA_works/works/lipeng/script/hallmark/mh.all.v2023.1.Mm.symbols.gmt",
    },
    "rat": {
        "go": "",
        "kegg": "/data/database/GSEA_gmt/rat/gene_kegg.gmt",
        "Hallmakr": "",
    },
}
