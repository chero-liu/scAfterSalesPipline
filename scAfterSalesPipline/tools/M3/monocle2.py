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
    --batch  F \\
    --subsampleid {self.subsampleid} \\
    --subgroup {self.subgroup} \\
    --subnew_celltype {self.subnew_celltype}

# vismethod: pseudotime,heatmap,treeplot,ridgeplot,module,trajectoryplot,expressplot_line,expressplot
    
Rscript  {self.script}/visualize_pseudotime.R  \\
    --input {self.outdir}/pseudotime_results.rds \\
    --genelist  ordering  \\
    --vismethod pseudotime \\
    --groupby clusters \\
    --branchpoint 1 \\
    --root 1 \\
    --CORES 5 \\
    --show_branch  {self.batch}  \\
    --anno {self.geneanno} \\
    --output {self.outdir}/visualize_pseudotime_results

    
# Rscript {self.script}/GeneSwitches.R \\
#     -i {self.input} \\
#     -m {self.outdir}/pseudotime_results.rds \\
#     -r {self.species} \\
#     -b 0.05 \\
#     --showtf T \\
#     -n 50 \\
#     -o /gpfs/oe-scrna/further_analysis/scRNA/10x/DZOE2024092169/result/M3/Monocle2/Microglia/visualize_pseudotime_results_root3_mergeState1_4/geneswitches


"""

        self.prefix = f"""{self.type}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
