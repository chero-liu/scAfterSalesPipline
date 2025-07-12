import os
import sys

from scAfterSalesPipline.tools.utils import ModuleFun


class Manualanno(ModuleFun):
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
        groupby: str,
    ):
        super().__init__(data, input, analysis, type)
        self.reduct = reduct
        self.celltype = celltype
        self.useColname = useColname
        self.extrabarcode = extrabarcode
        self.groupby = groupby

    def init_param(self):
        if self.reduct == None:
            self.reduct = "umap"
        if self.celltype == None:
            sys.exit("Error: celltype are required.")
        if self.groupby == None:
            self.groupby = "new_celltype"
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

Rscript  {self.script}/sctool \\
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

Rscript {self.script}/sctool \\
    --input  {self.outdir}/seurat.h5seurat \\
    --informat h5seurat \\
    --output {self.outdir} \\
    --assay RNA \\
    --dataslot data \\
    summarize \\
    --reduct {self.reduct} \\
    --palette customecol2 \\
    --groupby {self.groupby} \\
    --facetby sampleid,group \\
    --pointsize 0.5 \\
    --dosummary T \\
    --dims 2

Rscript  {self.script}/sctool \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/{self.groupby}  \\
    --outformat h5seurat   \\
    --ncores 10 \\
    --assay RNA  \\
    --dataslot counts,data,scale.data   \\
    --update T   \\
    findallmarkers \\
    --pct_fold 2 \\
    --topn_marker 10 \\
    --avg_log2FC 1.5 \\
    --pvalue 0.05 \\
    --FDR 0.1 \\
    --strict FALSE \\
    --test presto \\
    --cluster_name {self.groupby}


Rscript {self.script}/scVis \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/{self.groupby}  \\
    --ncores 10 \\
    --assay RNA \\
    --slot data,scale.data \\
    heatmap \\
    --markers {os.path.dirname(self.outdir)}/Markers/{self.groupby}/top10_markers_for_each_cluster.xls \\
    --topby gene_diff \\
    --topn 10 \\
    --groupby {self.groupby}  \\
    --group_colors customecol2 \\
    --sample_ratio 0.8 \\
    --style seurat 

Rscript {self.script}/sctool \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/{self.groupby}  \\
    --ncores 10 \\
    --assay RNA \\
    --dataslot data \\
    visualize \\
    --markers {os.path.dirname(self.outdir)}/Markers/{self.groupby}/top10_markers_for_each_cluster.xls \\
    --groupby {self.groupby} \\
    --reduct {self.reduct} \\
    --topn  10  \\
    --topby avg_log2FC \\
    --vismethod vlnplot,featureplot \\
    --vcolors customecol2 \\
    --ccolors spectral \\
    --dodge F

Rscript {self.script}/sctool annotation \\
  --genelist  {os.path.dirname(self.outdir)}/Markers/{self.groupby}/all_markers_for_each_cluster.xls \\
  --anno {self.geneanno}
Rscript {self.script}/sctool annotation \\
  --genelist {os.path.dirname(self.outdir)}/Markers/{self.groupby}/top10_markers_for_each_cluster.xls \\
  --anno {self.geneanno}
"""
        self.prefix = f"""{self.type}_{self.prefix}"""
        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
