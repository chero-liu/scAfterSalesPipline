import os
from pathlib import Path

from scAfterSalesPipline.tools.utils import ModuleFun


class Gsva(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        gmt: str,
        contrast: str,
        chunkby: str,
        downsample: str,
    ):
        super().__init__(data, input, analysis, type)
        self.gmt = gmt
        self.contrast = contrast
        self.chunkby = chunkby
        self.downsample = downsample

    def init_param(self):
        if self.chunkby == None:
            self.chunkby = 1000
        if self.downsample == None:
            self.downsample = 35000
        self.outdir = os.path.join(self.outdir, Path(self.gmt).stem)

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{self.environment}

Rscript /gpfs/oe-scrna/liuchenglong/RaD/scrnaAfterSalesPiplineScript/M2/gsva/GSVA_enrich.R \\
    --input {self.input} \\
    --gmt {self.gmt} \\
    --informat seurat   \\
    --OUTDIR  {self.outdir}  \\
    --method gsva \\
    --chunkby  {self.chunkby}  \\
    --kcdf Poisson \\
    --abs_rank  FALSE  \\
    --min_sz 2  \\
    --max_sz 10000  \\
    --parallel_sz 4 \\
    --mx_diff  TRUE \\
    --downsample {self.downsample} \\
    --prefix prefix \\
    --subnew_celltype {self.subnew_celltype} \\
    --subsampleid {self.subsampleid} \\
    --subgroup {self.subgroup}

for i in {" ".join(f'"{item}"' for item in self.contrast)}
do

y=$(echo "$i" | awk -F: '{{print $2 "_" $3}}')
z=$(echo "$i" | awk -F: '{{print $1}}')
Rscript /gpfs/oe-scrna/liuchenglong/RaD/scrnaAfterSalesPiplineScript/M2/gsva/GSVA_pathway_diffxp.R \\
    --input {self.outdir}/GSVA_enrichment_results.xls \\
    --seurat {self.input} \\
    --contrast $i \\
    --pval 0.05 \\
    --topn 10 \\
    --display_p TRUE \\
    --outdir {self.outdir}/$z/$y \\
    --prefix prefix \\
    --subnew_celltype {self.subnew_celltype} \\
    --subsampleid {self.subsampleid} \\
    --subgroup {self.subgroup}

done
"""

        self.prefix = f"""{self.type}_{Path(self.gmt).stem}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
