from scAfterSalesPipline.templates.path import *
from scAfterSalesPipline.tools.utils import Module1Fun


class FindMarker(Module1Fun):
    def __init__(
        self,
        data: dict,
        analysis: str,
        input: str,
        type: str,

        avg_log2FC: float,
        pvalue: float,
        useColname: str,
        reduct: str,
        Manualanno: str,
    ):
        super().__init__(data, input, analysis, type)

        self.avg_log2FC = avg_log2FC
        self.pvalue = pvalue
        self.useColname = useColname
        self.reduct = reduct
        self.Manualanno = Manualanno

    def init_param(self):
        if self.avg_log2FC == None:
            self.avg_log2FC = 1.5
        if self.pvalue == None:
            self.pvalue = 0.05
        if self.useColname is None and f"{self.Manualanno}" in f"{self.input}":
            self.useColname = "new_celltype"
        elif self.useColname == None:
            self.useColname = "clusters"
        if self.reduct == None:
            self.reduct = "umap"

        self.outdir = os.path.join(self.outdir, self.type,self.analysis,self.useColname)

    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

file_path={os.path.dirname(self.input)}/done.check

while true; do
    if [ -f "$file_path" ]; then
        break
    else
        sleep 60
    fi
done

Rscript  {self.analysis_script} \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}  \\
    --outformat h5seurat   \\
    --ncores 10 \\
    --assay RNA  \\
    --dataslot counts,data,scale.data   \\
    --update T   \\
    findallmarkers \\
    --pct_fold 2 \\
    --topn_marker 10 \\
    --avg_log2FC {self.avg_log2FC} \\
    --pvalue {self.pvalue} \\
    --FDR 0.1 \\
    --strict FALSE \\
    --test presto \\
    --cluster_name {self.useColname} 


Rscript {self.visualize_script} \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}  \\
    --ncores 10 \\
    --assay RNA \\
    --slot data,scale.data \\
    heatmap \\
    -l {self.outdir}/top10_markers_for_each_cluster.xls \\
    -c gene_diff \\
    -n 10 \\
    -g {self.useColname}  \\
    --group_colors customecol2 \\
    --sample_ratio 0.8 \\
    --style seurat 

Rscript {self.analysis_script} \\
    --input {self.input}  \\
    --informat h5seurat  \\
    --output {self.outdir}  \\
    --ncores 10 \\
    --assay RNA \\
    --dataslot data \\
    visualize \\
    -l {self.outdir}/top10_markers_for_each_cluster.xls \\
    -g {self.useColname} \\
    --reduct {self.reduct} \\
    --topn  10  \\
    --topby avg_log2FC \\
    -m vlnplot,featureplot \\
    --vcolors customecol2 \\
    --ccolors spectral \\
    --dodge F

Rscript {self.analysis_script} annotation \\
  --genelist  {self.outdir}/all_markers_for_each_cluster.xls \\
  --anno {GENE_ANNO[self.data["species"]]}
Rscript {self.analysis_script} annotation \\
  --genelist {self.outdir}/top10_markers_for_each_cluster.xls \\
  --anno {GENE_ANNO[self.data["species"]]}

"""

        self.prefix = f"""{self.prefix}_{self.type}"""
        self.save_script(shell_script_content)

    def run(self):
        self.init_param()
        self.shell_script()
