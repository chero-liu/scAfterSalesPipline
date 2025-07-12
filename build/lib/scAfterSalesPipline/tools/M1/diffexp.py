import os
import sys

from scAfterSalesPipline.tools.utils import ModuleFun


class Diffexp(ModuleFun):
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
        fdr: float = None,
        subobj: str = None,
        cluster_useColname: str = None,
    ):
        super().__init__(
            data,
            input,
            analysis,
            type,
        )

        self.contrast = contrast
        self.useColname = useColname
        self.FC = FC
        self.pvalue = pvalue
        self.fdr = fdr
        self.subobj = subobj
        self.cluster_useColname = cluster_useColname

    def init_param(self):
        if self.contrast == None:
            sys.exit("Error: contrast are required.")
        if self.useColname == None:
            self.useColname = "group"

        if self.cluster_useColname == None:
            self.cluster_useColname = "new_celltype"
        if self.FC == None:
            self.FC = 1.5
        if self.pvalue == None:
            self.pvalue = 0.05
        if self.fdr == None:
            self.fdr = ""
            self.p_name = "pval"
        else:
            self.fdr = f""" --fdr {self.fdr} """
            self.p_name = "padj"

        if self.subobj != None:
            tmp = f""" --predicate "{self.cluster_useColname} %in% c()" """
            for i in self.subobj:
                tmp = tmp.replace(")", f",\\'{i}\\')")
            self.predicate = tmp.replace("(,", "(")
            self.predicate_ = self.predicate + " &"
            self.outdir = os.path.join(
                self.outdir, "_".join(self.subobj), self.analysis
            )
            self.prefix = f"""{"_".join(self.subobj)}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}_{self.prefix}"""
        else:
            self.predicate = ""
            self.predicate_ = '--predicate ""'
            self.outdir = os.path.join(self.outdir, self.type, self.analysis)
            self.prefix = f"""{self.type}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}_{self.prefix}"""

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

Rscript  {self.script}/sctool \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}  \\
    --ncores 10 \\
    --assay RNA  \\
    --dataslot counts,data  {self.predicate}  \\
    diffexp \\
    --contrast {self.useColname}:{self.contrast} \\
    --FC {self.FC} \\
    --pvalue {self.pvalue} {self.fdr} \\
    --test presto

Rscript {self.script}/scVis \\
    --input  {self.input} \\
    --informat h5seurat \\
    --output {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]} \\
    --ncores 10 \\
    --assay RNA \\
    --slot data,scale.data {self.predicate_[:self.predicate_.rfind('"')]+self.predicate_[self.predicate_.rfind('"')+1:]}  {self.useColname} %in% c(\\'{self.contrast.split(":")[0]}\\',\\'{self.contrast.split(":")[1]}\\')" \\
    diff_heatmap \\
    --diffGene {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_*-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}.xls \\
    --topn 20 \\
    --groupby {self.useColname} \\
    --group_colors customecol2 \\
    --sample_ratio 0.8

Rscript {self.script}/sctool annotation \\
  --genelist  {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_*-all_diffexp_genes.xls \\
  --anno {self.geneanno}
Rscript {self.script}/sctool annotation \\
  --genelist  {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_*-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}.xls \\
  --anno {self.geneanno}


mkdir -p {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/enrichment_sh

perl  {os.path.dirname(self.script)}/enrich/enrich_go_kegg.pl  \\
-infile  {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_*-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}.xls \\
-go_bg {self.go_bg} \\
-category  {self.category}  \\
-kegg_bg  {self.kegg_bg}  \\
-outdir  {self.outdir}  \\
-shelldir {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/enrichment_sh \\
-thread 4 \\
-queue big

Rscript {os.path.dirname(self.script)}/volcanoplot/volcano.r \\
    -i {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-all_diffexp_genes.xls \\
    -P {self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]} \\
    -p {self.pvalue} {self.fdr} \\
    -f {self.FC} \\
    -o {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/Diffvolcano/ \\
    --symbol_topn 20


python {os.path.dirname(self.script)}/ppi/ppi.py \\
    --input {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}_anno.xls \\
    --species {self.ppispecies} \\
    --prefix {self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-{self.p_name}-{self.pvalue}-FC-{self.FC} \\
    --outdir {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/Diffppi/
Rscript  {os.path.dirname(self.script)}/ppi/ppi_circle.r  \\
   --input {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/Diffppi/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}.ppi_network.tsv \\
   --diff {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/{self.useColname}_{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}-diff-{self.p_name}-{self.pvalue}-FC-{self.FC}_anno.xls  \\
   --outpath {self.outdir}/{self.contrast.split(":")[0]}-vs-{self.contrast.split(":")[1]}/Diffppi/
"""

        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
