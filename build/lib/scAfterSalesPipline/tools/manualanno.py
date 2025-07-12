import os
import sys

from scAfterSalesPipline.tools.utils import Module1Fun


class Manualanno(Module1Fun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        reduct: str,
        celltype: str,
        useColname: str,
        extrabarcode: str,
    ):
        super().__init__(data, input, analysis, type)
        self.reduct = reduct
        self.celltype = celltype
        self.useColname = useColname
        self.extrabarcode = extrabarcode

    def init_param(self):
        if self.reduct == None:
            self.reduct = "umap"
        if self.celltype == None:
            sys.exit("Error: celltype are required.")
        if self.celltype.endswith(".tsv"):
            self.barcode = "F"
        elif self.celltype.endswith("csv"):
            self.barcode = "T"
        else:
            sys.exit("Error: celltype file format is not in [csv,tsv]")
        if self.useColname == None:
            self.useColname = "clusters"
        if self.extrabarcode != None:
            self.extrabarcode = f""" --extrabarcode "{self.extrabarcode}" """
        else:
            self.extrabarcode = ""

        self.outdir = os.path.join(self.outdir, self.type, self.analysis)

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

Rscript  {self.analysis_script} \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}  \\
    --outformat h5seurat   \\
    --ncores 10 \\
    --assay RNA  \\
    --dataslot counts,data,scale.data   \\
    --update F   \\
    changecelltype \\
    --celltype {self.celltype} \\
    --barcode {self.barcode} \\
    --cluster {self.useColname} \\
    --palette customecol2 \\
    --reduct {self.reduct} {self.extrabarcode}

touch {self.outdir}/done.check

Rscript {self.analysis_script} \\
    --input  {self.outdir}/seurat.h5seurat \\
    --informat h5seurat \\
    --output {self.outdir} \\
    --assay RNA \\
    --dataslot data \\
    summarize \\
    --reduct {self.reduct} \\
    --palette customecol2 \\
    --groupby {self.useColname} \\
    --facetby sampleid,group \\
    --pointsize 0.5 \\
    --dosummary T \\
    --dims 2

Rscript {self.visualize_script} \\
  --input {self.outdir}/seurat.h5seurat \\
  --informat h5seurat \\
  --output {self.outdir}/Correlation \\
  -t 6 \\
  --assay RNA \\
  --slot data \\
  --reduct {self.reduct} \\
  coefficient \\
  -g {self.useColname}

"""
        self.prefix = f"""{self.prefix}_{self.type}"""
        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
