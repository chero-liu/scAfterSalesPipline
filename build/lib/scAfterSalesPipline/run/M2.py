import os
import sys
import unittest

from scAfterSalesPipline.tools.M2.gsva import Gsva
from scAfterSalesPipline.tools.M2.addmodule import Addmodule
from scAfterSalesPipline.tools.utils import (
    Step,
    check_none,
    get_analysis,
    get_input,
    read_yaml,
    s_common,
)


class Module2(Step):
    def __init__(self, args):
        super().__init__(args)
        self.cwd = os.getcwd()
        self.input = (
            args.outdir + "/config/cfgM2.yaml"
            if args.input != "./config/cfgM2.yaml"
            else args.input
        )
        self.data = read_yaml(self.input)
        self.analysis_list = get_analysis(self.data["analysis"])
        self.module = "module"
        self.programID = "programID"
        self.species = "species"
        self.param = "param"
        self.h5seurat = "input"
        self.type = "type"
        self.prefix = "prefix"
        self.subnew_celltype = "subnew_celltype"
        self.subsampleid = "subsampleid"
        self.subgroup = "subgroup"
        self.Gsva = "GSVA"
        self.contrast = "contrast"
        self.gmt = "gmt"
        self.chunkby = "chunkby"
        self.downsample = "downsample"
        self.Addmodule = "Addmodule"
        self.extraGene = "extraGene"
        self.pvalue = "pvalue"
        self.scoredata = "scoredata"
        self.groupby = "groupby"

    def gsva(
        self,
        input=None,
        gmt=None,
        contrast=None,
    ):
        with Gsva(
            data=self.data,
            input=input,
            gmt=gmt,
            contrast=contrast,
            analysis=self.Gsva,
            type=self.data[self.param][self.type],
            chunkby=self.data[self.param][self.Gsva][self.chunkby],
            downsample=self.data[self.param][self.Gsva][self.downsample],
        ) as runner:
            runner.run()

    def addmodule(
        self,
        input=None,
    ):
        with Addmodule(
            data=self.data,
            input=input,
            extraGene=self.data[self.param][self.Addmodule][self.extraGene],
            pvalue=self.data[self.param][self.Addmodule][self.pvalue],
            scoredata=self.data[self.param][self.Addmodule][self.scoredata],
            groupby=self.data[self.param][self.Addmodule][self.groupby],
            analysis=self.Addmodule,
            type=self.data[self.param][self.type],
        ) as runner:
            runner.run()

    def run(self):
        if check_none(
            self.data[self.programID],
            self.data[self.species],
            self.data[self.param][self.type],
        ):
            sys.exit("Error: ProgramID , species and type are required.")

        if self.Gsva in self.analysis_list:
            input = get_input(
                analysis=self.Gsva,
                type=self.data[self.param][self.type],
                module=self.data[self.module],
                cwd=self.cwd,
                input=self.data[self.param][self.h5seurat],
            )

            if self.data[self.param][self.Gsva][self.gmt] == None:
                sys.exit("Error: gmt is required.")
            else:
                gmts = self.data[self.param][self.Gsva][self.gmt]

            if self.data[self.param][self.Gsva][self.contrast] == None:
                sys.exit("Error: contrast is required.")

            for gmt in gmts:
                self.gsva(
                    input=input,
                    gmt=gmt,
                    contrast=self.data[self.param][self.Gsva][self.contrast],
                )

        if self.Addmodule in self.analysis_list:
            input = get_input(
                analysis=self.Addmodule,
                type=self.data[self.param][self.type],
                module=self.data[self.module],
                cwd=self.cwd,
                input=self.data[self.param][self.h5seurat],
            )

            if self.data[self.param][self.Addmodule][self.extraGene] == None:
                sys.exit("Error: extraGene is required.")

            self.addmodule(
                input=input,
            )


def M2(args):
    with Module2(args) as runner:
        runner.run()


def get_opts_M2(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="./config/cfgM2.yaml",
        help="path of cfgM2.yaml",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
