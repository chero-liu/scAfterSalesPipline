import os
from pathlib import Path

from scAfterSalesPipline.templates.path import SCRIPT
from scAfterSalesPipline.tools.utils import Module2Fun


class Gsva(Module2Fun):
    def __init__(
        self,
        input: str,
        analysis: str,
        prefix: str,
        species: str,
        programID: str,
        type: str,
        outdir: str,
        subnew_celltype: str,
        subsampleid: str,
        subgroup: str,
        gmt: str,
        contrast: str,
        chunkby: int,
        downsample: int,
    ):
        self.cwd = os.getcwd()
        self.input = input
        self.analysis = analysis
        self.prefix = prefix
        self.species = species
        self.programID = programID
        self.type = type
        self.outdir = os.path.join(outdir, self.analysis, self.type)
        self.subnew_celltype = subnew_celltype
        self.subsampleid = subsampleid
        self.subgroup = subgroup

        self.gmt = gmt
        self.contrast = contrast
        self.chunkby = chunkby
        self.downsample = downsample

    def init_param(self):
        if self.chunkby == None:
            self.chunkby = 1000
        if self.downsample == None:
            self.downsample = 35000
        if self.subnew_celltype == None:
            self.subnew_celltype = "all"
        else:
            self.outdir = os.path.join(self.outdir, self.subnew_celltype)
            self.type = self.type + "_" + self.subnew_celltype
        if self.subsampleid == None:
            self.subsampleid = "all"
        else:
            self.outdir = os.path.join(self.outdir, self.subsampleid)
            self.type = self.type + "_" + self.subsampleid
        if self.subgroup == None:
            self.subgroup = "all"
        else:
            self.outdir = os.path.join(self.outdir, self.subgroup)
            self.type = self.type + "_" + self.subgroup
        self.outdir = os.path.join(self.outdir, Path(self.gmt).stem)

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{SCRIPT[self.analysis][2]}

Rscript {SCRIPT[self.analysis][0]} \\
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
Rscript {SCRIPT[self.analysis][1]} \\
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

        self.prefix = f"""{self.prefix}_{Path(self.gmt).stem}_{self.type}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
