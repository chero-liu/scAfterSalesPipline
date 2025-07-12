import os
from pathlib import Path

from scAfterSalesPipline.tools.utils import ModuleFun


class inferCNV(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        useColname: str,
        ref: str,
        clusting2use: str,
        refexp: str,
        downsample: str,
    ):
        super().__init__(data, input, analysis, type)
        self.useColname = useColname
        self.ref = ref
        self.clusting2use = clusting2use
        self.refexp = refexp
        self.downsample = downsample

    def init_param(self):
        if self.useColname == None:
            self.useColname = "new_celltype"
        if self.ref == None:
            self.ref = "ref"
        if self.clusting2use == None:
            self.clusting2use = "ward.D2"
        if self.refexp == None:
            self.refexp = ""
        else:
            self.refexp = f" --refexp {self.refexp}"
        if self.downsample == None:
            self.downsample = ""
        else:
            self.downsample = f" --down {self.downsample}"

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{self.environment}

Rscript {self.script}/infercnv.R \\
    --input {self.input} \\
    --informat seurat \\
    --celltype {self.useColname} \\
    --refgroup {self.ref} \\
    --gtf {self.gtf} \\
    --clusting2use  {self.clusting2use} {self.refexp} {self.downsample} \\
    --output {self.outdir}


Rscript {self.script}/infercnv_vis.R \\
    --input {self.input} \\
    --informat  seurat \\
    --output {self.outdir}/vis \\
    --infercnv_output_path {self.outdir} \\
    --groupby cnv_group,new_celltype,clusters,sampleid,group  \\
    --vismethod all \\
    --reduct umap \\
    --splitby sampleid,group \\
    --hmm F

"""

        self.prefix = f"""{self.type}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
