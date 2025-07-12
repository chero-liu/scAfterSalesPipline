import os
import unittest

from scAfterSalesPipline.templates.path import *
from scAfterSalesPipline.tools.utils import Step, s_common


class V4tov3(Step):
    def __init__(self, args):
        super().__init__(args)
        self.input = args.input

    def run(self):
        for input in self.input.split(","):
            os.system(
                f"""
                {ENVIRONMENT["v4tov3"]}
                Rscript {SCRIPT["v4tov3"]}/"v4tov3.r" \\
                    -i {input} \\
                    -f h5seurat \\
                    -o {os.path.dirname(input)}
                {ENVIRONMENT["GSVA"]}
                Rscript {SCRIPT["v4tov3"]}/"v4tov3.r" \\
                    -o {os.path.dirname(input)}
                """
            )


def v4tov3(args):
    with V4tov3(args) as runner:
        runner.run()


def get_opts_v4tov3(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="./",
        help="path of clean",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
