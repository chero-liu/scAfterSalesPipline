import os
from pathlib import Path

from scAfterSalesPipline.templates.path import SCRIPT
from scAfterSalesPipline.tools.utils import Module3Fun


class Monocle(Module3Fun):
    def __init__(
        self,
        input: str,
        analysis: str,
        prefix: str,
        species: str,
        programID: str,
        type: str,
        outdir: str,
        subobj: str,
        genelist: str,
    ):
        self.cwd = os.getcwd()
        self.input = input
        self.analysis = analysis
        self.prefix = prefix
        self.species = species
        self.programID = programID
        self.type = type
        self.outdir = outdir
        self.subobj = subobj
        self.genelist = genelist

    def init_param(self):
        pass

    def shell_script(self):
        shell_script_contenr = f"""
module purge
source /home/lipeng/miniconda3/bin/activate Scenic

# mkdir -p /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/
# cd /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/

# Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/scenic.R \
# -i /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M1/Macrophages/Clustering/data_ob_v3.rds \
# -f seurat \
# -d /data/database/SCENIC/mouse \
# -s mgi \
# --coexMethod top10perTarget \
# -o /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/

Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/RunRAS-RSS.R \
-i /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M1/Macrophages/Clustering/data_ob_v3.rds \
-v /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/int/3.4_regulonAUC.Rds \
-f rds \
-t 3 \
-c clusters \
-o /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/byclusters/RAS_RSS \
-s 0

Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/RunCSI.R \
-i /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M1/Macrophages/Clustering/data_ob_v3.rds \
-v /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/int/3.4_regulonAUC.Rds \
-f rds \
-c clusters \
-n 4 \
-o /gpfs/oe-scrna/liuchenglong/project/scRNA/DZOE2024011722/result/M6/scenic/Macrophages/byclusters/CSI \

"""

        self.prefix = f"""{self.prefix}_{self.type}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
