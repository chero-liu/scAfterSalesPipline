import os
from pathlib import Path
from scAfterSalesPipline.tools.utils import ModuleFun


class Markerplot(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        extraGene: str,
        groupby: str,
        reduct: str,
        vismethod: str,
    ):
        super().__init__(data, input, analysis, type)
        self.extraGene = extraGene
        self.groupby = groupby
        self.reduct = reduct
        self.vismethod = vismethod

    def init_param(self):
        if self.groupby == None:
            self.groupby = "clusters"
        if self.reduct == None:
            self.reduct = "umap"
        if self.vismethod == None:
            self.vismethod = "vlnplot,dotplot,featureplot"

        self.outdir = os.path.join(
            self.outdir,
            self.type,
            self.analysis,
            Path(self.extraGene).stem,
            self.groupby,
        )

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

Rscript  {self.script}/sctool \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}  \\
    --outformat h5seurat   \\
    --ncores 10 \\
    --assay RNA  \\
    --dataslot data \\
    visualize   \\
    --extraGene {self.extraGene} \\
    --groupby {self.groupby} \\
    --reduct {self.reduct} \\
    --vismethod {self.vismethod} \\
    --vcolors customecol2 \\
    --ccolors spectral \\
    --pointsize 0 \\
    --dodge F \\
    --dotsplit T

"""

        self.prefix = (
            f"""{self.type}_{Path(self.extraGene).stem}_{self.groupby}_{self.prefix}"""
        )
        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
