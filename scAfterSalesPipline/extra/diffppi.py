import os
import sys

from scAfterSalesPipline.templates.path import PPI_SPECIES
from scAfterSalesPipline.tools.utils import Module1Fun


class Diffppi(Module1Fun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        contrast: str = None,
        useColname: str = None,
        FC: float = None,
        pvalue: float = None,
        subobj: str = None,
    ):
        super().__init__(data, input, analysis, type)
        self.input = input
        self.contrast = contrast
        self.useColname = useColname
        self.FC = FC
        self.pvalue = pvalue
        self.subobj = subobj

    def init_param(self):
        if self.contrast == None:
            sys.exit("Error: contrast are required.")
        if self.useColname == None:
            self.useColname = "group"
        if self.FC == None:
            self.FC = 1.5
        if self.pvalue == None:
            self.pvalue = 0.05
        if self.subobj != None:
            self.input = os.path.join(self.outdir, "_".join(self.subobj), "Diffexp")
            self.outdir = os.path.join(
                self.outdir, "_".join(self.subobj), self.analysis
            )
            self.prefix = f"""{self.prefix}_{"_".join(self.subobj)}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}"""
        else:
            self.input = os.path.join(self.outdir, self.type, "Diffexp")
            self.outdir = os.path.join(self.outdir, self.type, self.analysis)
            self.prefix = f"""{self.prefix}_{self.type}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}"""

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

file_path={self.input}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/done.check

while true; do
    if [ -f "$file_path" ]; then
        break
    else
        sleep 60
    fi
done

python {self.analysis_script} \\
    --inputfile {self.input}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-*-{self.pvalue}-FC-{self.FC}_anno.xls \\
    --species {PPI_SPECIES[self.species]} \\
    --prefix {self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-*-{self.pvalue}-FC-{self.FC} \\
    --outputdir {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]} \\
    --noft 25
Rscript  {self.visualize_script}  \\
   --input {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-*-{self.pvalue}-FC-{self.FC}.ppi_network.tsv \\
   --diff {self.input}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-*-{self.pvalue}-FC-{self.FC}_anno.xls  \\
   --output {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}
"""

        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
