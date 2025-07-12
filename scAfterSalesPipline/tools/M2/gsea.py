import os
from pathlib import Path

from scAfterSalesPipline.tools.utils import ModuleFun


class Gsea(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        gmt: str,
        contrast: str,
    ):
        super().__init__(data, input, analysis, type)
        self.gmt = gmt
        self.contrast = contrast

    def init_param(self):
        self.outdir = os.path.join(self.outdir, Path(self.gmt).stem)

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{self.environment}
 
# export PYTHONPATH=/gpfs/oe-scrna/guopengyu/GSEA/oebio:$PYTHONPATH

for i in {" ".join(f'"{item}"' for item in self.contrast)}
do

y=$(echo "$i" | awk -F: '{{print $2 "_" $3}}')
z=$(echo "$i" | awk -F: '{{print $1}}')

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scAfterSalesPipline/scAfterSalesPipline/script/M2/gsea/GSEA_sc.R \\
    --input {self.input}  \\
    --contrast $i  \\
    --minsize 15 \\
    --maxsize  500  \\
    --gmt {self.gmt}  \\
    --outdir {self.outdir}/$z/$y  \\
    --subnew_celltype {self.subnew_celltype} \\
    --subsampleid {self.subsampleid} \\
    --subgroup {self.subgroup}  \\
    --subcluster {self.subcluster}

done
"""

        self.prefix = f"""{self.type}_{Path(self.gmt).stem}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
