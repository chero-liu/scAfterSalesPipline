import os
import sys
import unittest

from scAfterSalesPipline.tools.M1.cluster import Cluster
from scAfterSalesPipline.tools.M1.enrich import Enrich
from scAfterSalesPipline.tools.M1.diffexp import Diffexp
from scAfterSalesPipline.tools.M1.manualanno import Manualanno
from scAfterSalesPipline.tools.M1.markerplot import Markerplot
from scAfterSalesPipline.tools.utils import (
    Step,
    check_none,
    get_analysis,
    get_input,
    read_yaml,
    s_common,
)


class Module1(Step):
    def __init__(self, args):
        super().__init__(args)
        self.cwd = os.getcwd()
        self.input = (
            args.outdir + "/config/cfgM1.yaml"
            if args.input != "./config/cfgM1.yaml"
            else args.input
        )
        self.data = read_yaml(self.input)
        self.analysis_list = get_analysis(self.data["analysis"])
        self.module = self.data["module"]
        self.programID = self.data["programID"]
        self.species = self.data["species"]

        self.h5seurat = self.data["param"]["input"]
        self.type = self.data["param"]["type"]
        self.prefix = self.data["param"]["prefix"]
        self.Clustering = "Clustering"
        self.subobj = self.data["param"]["Clustering"]["subobj"]
        self.useColname = self.data["param"]["Clustering"]["useColname"]
        self.reduct1 = self.data["param"]["Clustering"]["reduct1"]
        self.reduct2 = self.data["param"]["Clustering"]["reduct2"]
        self.batchid = self.data["param"]["Clustering"]["batchid"]
        self.resolution = self.data["param"]["Clustering"]["resolution"]
        self.Manualanno = "Manualanno"
        self.celltype = self.data["param"]["Manualanno"]["celltype"]
        self.extrabarcode = self.data["param"]["Manualanno"]["extrabarcode"]
        self.orderby = self.data["param"]["Manualanno"]["orderby"]
        self.Diffexp = "Diffexp"
        self.contrast = self.data["param"]["Diffexp"]["contrast"]
        self.FC = self.data["param"]["Diffexp"]["FC"]
        self.pvalue = self.data["param"]["Diffexp"]["pvalue"]
        self.fdr = self.data["param"]["Diffexp"]["fdr"]
        self.Enrich = "Enrich"
        self.PPI = "PPI"
        self.splitby = self.data["param"]["Enrich"]["splitby"]
        self.sortby = self.data["param"]["Enrich"]["sortby"]
        self.top = self.data["param"]["Enrich"]["top"]
        self.last = self.data["param"]["Enrich"]["last"]
        self.Markerplot = "Markerplot"
        self.extraGene = self.data["param"][self.Markerplot]["extraGene"]
        self.groupby = self.data["param"][self.Markerplot]["groupby"]
        self.vismethod = self.data["param"][self.Markerplot]["vismethod"]

    def cluster(self, input=None, subobj=None, type=None):
        with Cluster(
            data=self.data,
            input=input,
            analysis=self.Clustering,
            type=type,
            subobj=subobj,
            batchid=self.batchid,
            resolution=self.resolution,
            reduct1=self.reduct1,
            reduct2=self.reduct2,
            useColname=self.data["param"]["Clustering"]["useColname"],
        ) as runner:
            runner.run()

    def manualanno(self, input=None):
        with Manualanno(
            data=self.data,
            input=input,
            type=self.type,
            analysis=self.Manualanno,
            reduct=self.reduct2,
            celltype=self.celltype,
            extrabarcode=self.extrabarcode,
            groupby=None,
            orderby=self.orderby,
        ) as runner:
            runner.run()

    def markerplot(self, input=None, type=None):
        with Markerplot(
            data=self.data,
            input=input,
            analysis=self.Markerplot,
            type=type,
            extraGene=self.extraGene,
            groupby=self.groupby,
            reduct=self.reduct2,
            vismethod=self.vismethod,
        ) as runner:
            runner.run()

    def diffexp(self, input=None, contrast=None, subobj=None, type=None):
        with Diffexp(
            data=self.data,
            input=input,
            contrast=contrast,
            analysis=self.Diffexp,
            useColname=self.data["param"]["Diffexp"]["useColname"],
            cluster_useColname=self.data["param"]["Clustering"]["useColname"],
            FC=self.FC,
            pvalue=self.pvalue,
            fdr=self.fdr,
            subobj=subobj,
            type=type,
        ) as runner:
            return runner.run()

    def enrich(self, input=None, type=None):
        with Enrich(
            data=self.data,
            input=input,
            analysis=self.Enrich,
            type=type,
            splitby=self.splitby,
            sortby=self.sortby,
            top=self.top,
            last=self.last,
        ) as runner:
            runner.run()

    def run(self):
        # check required parameters
        if check_none(
            self.programID,
            self.species,
            self.type,
        ):
            sys.exit("Error: ProgramID , species and type are required.")

        # Generate shell scripts
        if self.Clustering in self.analysis_list:
            if isinstance(self.subobj, list):
                for i in range(len(self.subobj)):
                    subobj = self.subobj[i]
                    input = get_input(
                        self.Clustering,
                        self.module,
                        self.cwd,
                        self.h5seurat,
                        subobj,
                        type=self.type,
                        useColname=self.data["param"]["Clustering"]["useColname"],
                    )
                    self.cluster(input, subobj, type=self.type)

            elif self.subobj == None:
                input = get_input(
                    self.Clustering,
                    self.module,
                    self.cwd,
                    self.h5seurat,
                    self.subobj,
                    type=self.type,
                    useColname=self.data["param"]["Clustering"]["useColname"],
                )
                self.cluster(input, self.subobj, type=self.type)
            else:
                print("Error: subobj must be a list.")

        if self.Manualanno in self.analysis_list:
            input = get_input(
                self.Manualanno,
                self.module,
                self.cwd,
                self.h5seurat,
                self.subobj,
                type=self.type,
            )
            self.manualanno(input=input)

        if self.Markerplot in self.analysis_list:
            input = get_input(
                self.Markerplot,
                self.module,
                self.cwd,
                self.h5seurat,
                self.subobj,
                type=self.type,
            )
            self.markerplot(input, type=self.type)

        if self.Diffexp in self.analysis_list:
            if isinstance(self.contrast, list) and isinstance(self.subobj, list):
                for contrast in self.contrast:
                    for subobj in self.subobj:
                        input = get_input(
                            self.Diffexp,
                            self.module,
                            self.cwd,
                            self.h5seurat,
                            "_".join(subobj),
                            type=self.type,
                            useColname=self.data["param"]["Clustering"]["useColname"],
                        )
                        self.diffexp(input, contrast, subobj, type=self.type)
            elif isinstance(self.contrast, list) and self.subobj == None:
                for contrast in self.contrast:
                    input = get_input(
                        self.Diffexp,
                        self.module,
                        self.cwd,
                        self.h5seurat,
                        self.subobj,
                        type=self.type,
                    )
                    self.diffexp(input, contrast, self.subobj, type=self.type)
            else:
                print("Error: contrast must be a list.")

        if self.Enrich in self.analysis_list:
            input = get_input(
                self.Enrich,
                self.module,
                self.cwd,
                self.h5seurat,
                self.subobj,
                type=self.type,
            )
            self.enrich(input, type=self.type)


def M1(args):
    with Module1(args) as runner:
        runner.run()


def get_opts_M1(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="./config/cfgM1.yaml",
        help="path of cfgM1.yaml",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
