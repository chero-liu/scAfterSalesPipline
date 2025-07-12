import os
from pathlib import Path

from scAfterSalesPipline.tools.utils import ModuleFun


class Addmodule(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        extraGene: str,
        pvalue: str,
        scoredata: str,
        groupby: str,
    ):
        super().__init__(data, input, analysis, type)
        self.extraGene = extraGene
        self.pvalue = pvalue
        self.scoredata = scoredata
        self.groupby = groupby

    def init_param(self):
        if self.groupby == None:
            self.groupby = "clusters"
        if self.pvalue == None:
            self.pvalue = ""
        else:
            self.pvalue = f" --pvalue {self.pvalue}"

        if self.scoredata == None:
            self.scoredata = ""
        else:
            self.scoredata = f" --scoredata {self.scoredata}"
        self.outdir = os.path.join(self.outdir, Path(self.extraGene).stem)

    def shell_script(self):
        shell_script_contenr = f"""
#!/bin/bash
{self.environment}

Rscript   {self.script}/visualize_markers.R  \\
    --RDS {self.input} \\
    --extraGene {self.extraGene} \\
    --output {self.outdir} \\
    --vismethod geneset  \\
    --pointsize 0.5 \\
    --assay RNA \\
    --reduct umap  \\
    --groupby {self.groupby}  {self.pvalue}  {self.scoredata}

"""

        self.prefix = f"""{self.type}_{Path(self.extraGene).stem}_{self.prefix}"""
        self.save_script(shell_script_contenr)

    def run(self):
        self.init_param()
        self.shell_script()
