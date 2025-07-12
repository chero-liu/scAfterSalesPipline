import os
import sys

from scAfterSalesPipline.tools.utils import ModuleFun


class Cluster(ModuleFun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,
        subobj: str,
        batchid: str,
        resolution: str,
        reduct1: str,
        reduct2: str,
        useColname: str,
    ):
        super().__init__(data, input, analysis, type)
        self.subobj = subobj
        self.batchid = batchid
        self.resolution = resolution
        self.reduct1 = reduct1
        self.reduct2 = reduct2
        self.useColname = useColname

    def init_param(self):
        if self.resolution == None:
            self.resolution = 0.4

        if self.reduct2 == None:
            self.reduct2 = "umap"

        if self.batchid == None:
            self.batchid = "batchid"

        if self.useColname == None:
            self.useColname = "new_celltype"

        if self.subobj != None:
            tmp = f""" --predicate "{self.useColname} %in% c()" """
            for i in self.subobj:
                tmp = tmp.replace(")", f",\\'{i}\\')")
            self.predicate = tmp.replace("(,", "(")
            self.annolevel = "single"
            self.outdir = os.path.join(
                self.outdir, "_".join(self.subobj), self.analysis
            )
        else:
            self.predicate = ""
            self.annolevel = "main"
            self.outdir = os.path.join(self.outdir, self.type, self.analysis)
        if self.reduct1 == None or self.reduct1 == "pca":
            self.reduct1 = "pca"
        elif self.reduct1 == "mnn":
            self.reduct1 = f"{self.reduct1}   --components 10 "
        elif self.reduct1 == "harmony":
            self.reduct1 = f"pca,{self.reduct1}   -t 20 -y 30 "
        else:
            sys.exit("reduct1 is not in [mnn,harmony,pca]")

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
    --update F  {self.predicate} \\
    bclust   \\
    --reduct1 {self.reduct1}  \\
    --batchid {self.batchid} \\
    --reduct2 {self.reduct2}   \\
    --clusteringuse snn  \\
    --resolution {self.resolution}   \\
    --rerun T   \\
    --pointsize  0.5   \\
    --palette customecol2

Rscript {self.script}/sctool \\
    --input  {self.outdir}/seurat.h5seurat \\
    --informat h5seurat \\
    --output {self.outdir} \\
    --assay RNA \\
    --dataslot data \\
    summarize \\
    --reduct {self.reduct2} \\
    --palette customecol2 \\
    --groupby clusters \\
    --facetby sampleid,group \\
    --pointsize 0.5 \\
    --dosummary T

Rscript {self.script}/scVis \\
    --input {self.outdir}/seurat.h5seurat \\
    --informat h5seurat \\
    --output {self.outdir}/Correlation \\
    -t 6 \\
    --assay RNA \\
    --slot data \\
    --reduct {self.reduct2} \\
    coefficient \\
    -g clusters

Rscript  {self.script}/sctool  \\
    --input {self.outdir}/seurat.h5seurat \\
    --informat h5seurat \\
    --output {os.path.dirname(self.outdir)}/Reference_celltype \\
    --outformat h5seurat \\
    --update T \\
    --assay RNA \\
    --dataslot counts \\
    celltyping \\
    --customref {self.cellanno} \\
    --annolevel {self.annolevel} \\
    --usecluster F \\
    --demethod classic \\
    --pointsize 0.3 \\
    --topn 25 \\
    --reduct {self.reduct2} \\
    --species {self.species}

Rscript  {self.script}/sctool \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/clusters/  \\
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
    --cluster_name clusters


Rscript {self.script}/scVis \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/clusters  \\
    --ncores 10 \\
    --assay RNA \\
    --slot data,scale.data \\
    heatmap \\
    -l {os.path.dirname(self.outdir)}/Markers/clusters/top10_markers_for_each_cluster.xls \\
    -c gene_diff \\
    -n 10 \\
    -g clusters  \\
    --group_colors customecol2 \\
    --sample_ratio 0.8 \\
    --style seurat 

Rscript {self.script}/sctool \\
    --input {self.outdir}/seurat.h5seurat  \\
    --informat h5seurat  \\
    --output {os.path.dirname(self.outdir)}/Markers/clusters  \\
    --ncores 10 \\
    --assay RNA \\
    --dataslot data \\
    visualize \\
    -l {os.path.dirname(self.outdir)}/Markers/clusters/top10_markers_for_each_cluster.xls \\
    -g clusters \\
    --reduct {self.reduct2} \\
    --topn  10  \\
    --topby avg_log2FC \\
    -m vlnplot,featureplot \\
    --vcolors customecol2 \\
    --ccolors spectral \\
    --dodge F

Rscript {self.script}/sctool annotation \\
  --genelist  {os.path.dirname(self.outdir)}/Markers/clusters/all_markers_for_each_cluster.xls \\
  --anno {self.geneanno}
Rscript {self.script}/sctool annotation \\
  --genelist {os.path.dirname(self.outdir)}/Markers/clusters/top10_markers_for_each_cluster.xls \\
  --anno {self.geneanno}

"""
        if self.subobj != None:
            self.prefix = f"""{"_".join(self.subobj)}_{self.prefix}"""
        else:
            self.prefix = f"""{self.type}_{self.prefix}"""
        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
