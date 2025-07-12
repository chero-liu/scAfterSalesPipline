import os
import sys
import unittest

from scAfterSalesPipline.tools.M3.monocle2 import Monocle2
from scAfterSalesPipline.tools.utils import (
    Step,
    check_none,
    get_analysis,
    get_input,
    read_yaml,
    s_common,
)


class Module3(Step):
    def __init__(self, args):
        super().__init__(args)
        self.cwd = os.getcwd()
        self.input = (
            args.outdir + "/config/cfgM3.yaml"
            if args.input != "./config/cfgM3.yaml"
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
        self.Monocle2 = "Monocle2"
        self.useColname = "useColname"
        self.downsample = "downsample"
        self.colorby = "colorby"
        self.resolution = "resolution"
        self.batch = "batch"

    def monocle2(
        self,
        input=None,
    ):
        with Monocle2(
            data=self.data,
            input=input,
            analysis=self.Monocle2,
            type=self.data[self.param][self.type],
            useColname=self.data[self.param][self.Monocle2][self.useColname],
            downsample=self.data[self.param][self.Monocle2][self.downsample],
            colorby=self.data[self.param][self.Monocle2][self.colorby],
            resolution=self.data[self.param][self.Monocle2][self.resolution],
            batch=self.data[self.param][self.Monocle2][self.batch],
        ) as runner:
            runner.run()

    def run(self):
        if check_none(
            self.data[self.programID],
            self.data[self.species],
            self.data[self.param][self.type],
        ):
            sys.exit("Error: ProgramID , species and type are required.")

        if self.Monocle2 in self.analysis_list:
            input = get_input(
                analysis=self.Monocle2,
                type=self.data[self.param][self.type],
                module=self.data[self.module],
                cwd=self.cwd,
                input=self.data[self.param][self.h5seurat],
            )

            self.monocle2(input=input)


def M3(args):
    with Module3(args) as runner:
        runner.run()


def get_opts_M3(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="./config/cfgM3.yaml",
        help="path of cfgM3.yaml",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
