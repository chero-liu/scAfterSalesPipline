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

#mgi  /data/database/SCENIC/mouse
output_dir=/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024102138/result/M6/Scenic/Like-macrophage
input=/gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024102138/result/M1/Like-macrophage/Manualanno/data_ob_v3.rds
group=group

mkdir -p $output_dir
cd $output_dir
Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M6/scenic/scenic.R \
    -i $input \
    -f seurat \
    -d /data/database/SCENIC/human/ \
    -s hgnc \
    --coexMethod top10perTarget \
    -o $output_dir

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M6/scenic/RunRAS-RSS.R \
    -i $input \
    -v int/3.4_regulonAUC.Rds \
    -f rds \
    -t 3 \
    -c $group \
    -o $output_dir/by$group \
    -s 0

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M6/scenic/RunCSI.R \
    -i $input \
    -v int/3.4_regulonAUC.Rds \
    -f rds \
    -c $group \
    -n 4 \
    -o $output_dir/by$group
"""

        self.prefix = f"""{self.prefix}_{self.type}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
