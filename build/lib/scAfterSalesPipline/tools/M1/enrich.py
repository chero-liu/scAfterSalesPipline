import os

from scAfterSalesPipline.tools.utils import ModuleFun


class Enrich(ModuleFun):
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
    ):
        super().__init__(data, input, analysis, type)
        self.input = input
        self.splitby = splitby
        self.sortby = sortby
        self.top = top
        self.last = last

    def init_param(self):
        if self.splitby == None:
            self.splitby = "'None'"
        if self.sortby == None:
            self.sortby = "'None'"
        if self.top == None:
            self.top = "'None'"
        if self.last == None:
            self.last = "'None'"

        self.outdir = os.path.join(self.outdir, self.type, self.analysis, self.prefix)
        self.prefix = f"""{self.type}_{self.prefix}"""

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

python {self.script}/splitByCol.py \\
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

for i in `ls {self.outdir}/${{x}}*.xls`
do
perl  {self.script}/enrich_go_kegg.pl  \\
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
