import os
import sys

from scAfterSalesPipline.templates.path import ENRICH_BG
from scAfterSalesPipline.tools.utils import Module1Fun


class Diffenrich(Module1Fun):
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
        go_bg: str = None,
        category: str = None,
        kegg_bg: str = None,
    ):
        super().__init__(data, input, analysis, type)
        self.input = input
        self.contrast = contrast
        self.go_bg = go_bg
        self.category = category
        self.kegg_bg = kegg_bg
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
        if self.go_bg == None:
            self.go_bg = ENRICH_BG[self.species]["go_bg"]
        if self.category == None:
            self.category = ENRICH_BG[self.species]["category"]
        if self.kegg_bg == None:
            self.kegg_bg = ENRICH_BG[self.species]["kegg_bg"]

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

mkdir -p {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/enrichment_sh

perl  {self.visualize_script}  \\
-infile  {self.input}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_*-diff-*-{self.pvalue}-FC-{self.FC}.xls \\
-go_bg {self.go_bg} \\
-category  {self.category}  \\
-kegg_bg  {self.kegg_bg}  \\
-outdir  {self.outdir}  \\
-shelldir {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/enrichment_sh \\
-thread 4 \\
-queue big
"""

        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
