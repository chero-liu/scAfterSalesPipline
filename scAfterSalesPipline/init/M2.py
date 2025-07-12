import os
import sys
import unittest

from ruamel.yaml import YAML
from scAfterSalesPipline.init.__init__ import M2DIRECTORIES
from scAfterSalesPipline.templates.path import CONFIG
from scAfterSalesPipline.tools.utils import Step, s_common


class Module2(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.prefix = args.prefix
        self.input = args.input
        self.outdir = args.outdir
        self.species = args.species
        self.type = args.type
        self.directory_list = M2DIRECTORIES

    def init_dirs(self):
        for directory in self.directory_list:
            dir_path = os.path.join(self.outdir, directory)
            try:
                os.makedirs(dir_path, exist_ok=True)
            except Exception as e:
                print(f"Failed to create directory {dir_path}: {e}", file=sys.stderr)
                sys.exit(1)

        print("Directory structure initialized successfully")

    def init_yaml(self):
        yaml = YAML()

        with open(
            CONFIG["M2"],
            "r",
        ) as file:
            data = yaml.load(file)

        if self.prefix:
            data["programID"] = self.prefix
        if self.species:
            data["species"] = self.species

        with open(os.path.join(self.outdir, "config", "cfgM2.yaml"), "w") as file:
            yaml.dump(data, file)

    def load_data(self):
        type_path = os.path.join(self.outdir, "result", "M2", self.type)
        if not os.path.exists(type_path):
            os.makedirs(os.path.exists(type_path))

        os.system(f"cp {self.input} {type_path}/{self.prefix}.h5seurat")

    def run(self):
        self.init_dirs()
        self.init_yaml()
        if self.input:
            self.load_data()


def M2(args):
    with Module2(args) as runner:
        runner.run()


def get_opts_M2(parser, sub_program=True):
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="path to the seurat.h5seurat dir",
    )

    parser.add_argument(
        "--type",
        type=str,
        help="type of the input seurat.h5seurat",
    )

    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
