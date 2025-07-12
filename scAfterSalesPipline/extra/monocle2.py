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
#!/bin/bash

module load OESingleCell/3.0.d
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
    --input {self.input}  \\
    --informat h5seurat \\
    --assay RNA \\
    --output {self.outdir} \\
    --ncores 8 \\
    --update FALSE \\
    monocle \\
    --design {self.design} \\
    --min.cell 0.01 \\
    --resolution 0.4 \\
    --colorby {self.colorby} \\
    --pointsize 1

Rscript   /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \\
    --input {self.outdir}/pseudotime_results.rds \\
    --output {self.outdir}/visualize_pseudotime_results \\
    --ncores 6 \\
    monocle_vis \\
    --genelist  {self.genelist}  \\
    --vismethod {self.vismethod} \\
    --groupby {self.groupby}  \\
    --show_branch  FALSE \\
    --anno  \\
    --root 1 \\
    --branchpoint 1
"""

        self.prefix = f"""{self.prefix}_{self.type}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
