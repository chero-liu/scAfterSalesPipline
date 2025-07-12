import os

from scAfterSalesPipline.tools.utils import ModuleFun


class Monocle2(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        useColname: str,
        resolution: str,
        colorby: str,
        downsample: str,
        batch: str,
    ):
        super().__init__(data, input, analysis, type)
        self.useColname = useColname
        self.resolution = resolution
        self.colorby = colorby
        self.downsample = downsample
        self.batch = batch

    def init_param(self):
        if self.useColname == None:
            self.useColname = "clusters"
        if self.resolution == None:
            self.resolution = 0.4
        if self.colorby == None:
            self.colorby = "clusters,sampleid,group"
        if self.downsample == None:
            self.downsample = 30000
        if self.batch == None:
            self.batch = "F"

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{self.environment}

Rscript {self.script}/monocle.R \\
    --INEXPRESS {self.input}  \\
    --INFORMAT seurat \\
    --assay RNA \\
    --design {self.useColname} \\
    --CORES 8 \\
    --min.cell 0.01 \\
    --resolution {self.resolution} \\
    --colorby {self.colorby} \\
    --pointsize 1 \\
    --output {self.outdir} \\
    --downsample {self.downsample} \\
    --batch  {self.batch} \\
    --subsampleid {self.subsampleid} \\
    --subgroup {self.subgroup} \\
    --subnew_celltype {self.subnew_celltype}


# Rscript  /public/scRNA_works/pipeline/scRNA-seq_further_analysis/visualize_pseudotime.R  \\
#     -i pseudotime_results.rds \\
#     -g  ordering  \\
#     -m all \\
#     -c clusters \\
#     -b 1 \\
#     -j 5 \\
#     -a /data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls \\
#     -o visualize_pseudotime_results

"""

        self.prefix = f"""{self.type}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
