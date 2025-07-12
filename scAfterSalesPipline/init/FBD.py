import os
import sys
import unittest
import re
from scAfterSalesPipline.tools.utils import Step, s_common
import nbformat


from datetime import datetime

date = datetime.today().strftime("%Y-%m-%d")


class FeedbackResultsDirectory(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.prefix = args.prefix
        self.input = args.input
        self.outdir = args.outdir
        self.species = args.species
        self.type = args.type
        self.date = date
        self.date_format = r"^\d{4}-\d{2}-\d{2}$"

    def check_directory_for_date_format(self):

        for item in os.listdir(self.outdir):
            item_path = os.path.join(self.outdir, item)

            if os.path.isdir(item_path) and re.match(self.date_format, item):
                return os.path.join(self.outdir, item, "result")

        return False

    def update_directory_date(self):
        self.check_result = self.check_directory_for_date_format()

        if self.check_result == False:
            return self.check_result
        else:
            pattern = r"^(.*)\d{4}-\d{2}-\d{2}$"
            for item in os.listdir(self.check_result):
                item_path = os.path.join(self.check_result, item)
                if os.path.isdir(item_path) and re.match(pattern, item):
                    updated_item = re.sub(r"\d{4}-\d{2}-\d{2}$", self.date, item)
                    return updated_item
                else:
                    return None

    def init_dirs(self):
        if self.input != None:
            self.date = self.input

        self.directory_list = [
            f"{self.date}/data",
            f"{self.date}/result",
            f"{self.date}/script",
        ]

        update_result = self.update_directory_date()

        for directory in self.directory_list:
            dir_path = os.path.join(self.outdir, directory)
            try:
                os.makedirs(dir_path, exist_ok=True)
            except Exception as e:
                print(f"Failed to create directory {dir_path}: {e}", file=sys.stderr)
                sys.exit(1)

        if update_result == False:
            print("No directory xxxx-xx-xx", file=sys.stderr)
        elif update_result == None:
            print("No directory yyyyyyyyyyxxxx-xx-xx", file=sys.stderr)
        else:
            os.makedirs(
                os.path.join(self.outdir, self.date, "result", update_result),
                exist_ok=True,
            )

            with open(
                os.path.join(
                    self.outdir, self.date, "result", update_result, "README.txt"
                ),
                "w",
            ) as f:
                f.write(
                    f"# {self.prefix}\n\n"
                    f"## {self.type}\n\n"
                    f"## {self.species}\n\n"
                    f"## {self.date}\n\n"
                )

            nb = nbformat.v4.new_notebook()
            code_cell = nbformat.v4.new_code_cell(
                "source('~/script/lclFunc.r')\nlibrary(Seurat)\nlibrary(ggplot2)\nlibrary(dplyr)\nlibrary(tidyr)\n\nimport numpy as np\nimport pandas as pd\nimport matplotlib.pyplot as plt\nimport seaborn as sns\nimport os\nimport sys\nimport re\nsys.path.append('~/script/lclFunc.py')"
            )
            nb.cells.append(code_cell)

            with open(
                os.path.join(self.outdir, self.date, "script", "test.ipynb"),
                "w",
                encoding="utf-8",
            ) as f:
                nbformat.write(nb, f)

    def run(self):
        self.init_dirs()


def FBD(args):
    with FeedbackResultsDirectory(args) as runner:
        runner.run()


def get_opts_FBD(parser, sub_program=True):
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
