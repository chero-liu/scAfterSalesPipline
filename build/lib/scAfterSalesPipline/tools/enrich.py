import os
import sys

from scAfterSalesPipline.templates.path import ENRICH_BG
from scAfterSalesPipline.tools.utils import Module1Fun


class Enrich(Module1Fun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        splitby: str,
        sortby: str,
        top: int,
        last: int,
        go_bg: str = None,
        category: str = None,
        kegg_bg: str = None,
    ):
        super().__init__(data, input, analysis, type)
        self.input = input
        self.splitby = splitby
        self.sortby = sortby
        self.top = top
        self.last = last
        self.go_bg = go_bg
        self.category = category
        self.kegg_bg = kegg_bg


    def init_param(self):
        if self.splitby == None:
            self.splitby = "'None'"
        if self.sortby == None:
            self.sortby = "'None'"
        if self.top == None:
            self.top = "'None'"
        if self.last == None:
            self.last = "'None'"
        if self.go_bg == None:
            self.go_bg = ENRICH_BG[self.species]["go_bg"]
        if self.category == None:
            self.category = ENRICH_BG[self.species]["category"]
        if self.kegg_bg == None:
            self.kegg_bg = ENRICH_BG[self.species]["kegg_bg"]

        self.outdir = os.path.join(self.outdir, self.type, self.analysis, self.prefix)
        self.prefix = f"""{self.prefix}_{self.type}"""

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

python {self.analysis_script} \\
    --input {self.input}  \\
    --outdir {self.outdir} \\
    --splitby {self.splitby} \\
    --sortby {self.sortby} \\
    --top {self.top} \\
    --last {self.last}

x={self.splitby}

if [ "$x" == "None" ]; then
  x=""
fi

for i in `ls {self.outdir}/${{x}}*.xls`;do
perl  {self.visualize_script}  \\
    -infile  $i \\
    -go_bg {self.go_bg} \\
    -category  {self.category}  \\
    -kegg_bg  {self.kegg_bg}  \\
    -outdir  {self.outdir}  \\
    -shelldir {self.outdir}/enrichment_sh \\
    -thread 4 \\
    -queue big
done
"""

        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
