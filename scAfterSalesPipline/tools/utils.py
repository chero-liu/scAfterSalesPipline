import abc
import importlib
import os
import sys
from pathlib import Path

import yaml
from scAfterSalesPipline.__init__ import ROOT_PATH
from scAfterSalesPipline.templates.path import *


def read_yaml(file_path):
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)
    return data


def check_none(*args):
    return any(arg is None for arg in args)


def get_analysis(dic, condition=lambda x: x != 0):
    result = [key for key, value in dic.items() if condition(value)]
    if not result:
        print("No analysis found, Please select at least one valid analysis")
        sys.exit(1)
    else:
        print(f"Analysis list: {result}")
    return result


def find_assay_init(assay):
    init_module = importlib.import_module(f"scAfterSalesPipline.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    file_path_dict = {
        "assay": f"{ROOT_PATH}/{assay}/{step}.py",
        "tools": f"{ROOT_PATH}/tools/{step}.py",
    }

    init_module = find_assay_init(assay)

    if os.path.exists(file_path_dict["assay"]):
        step_module = importlib.import_module(f"scAfterSalesPipline.{assay}.{step}")
    elif hasattr(init_module, "IMPORT_DICT") and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict["tools"]):
        step_module = importlib.import_module(f"scAfterSalesPipline.tools.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module


def s_common(parser):
    """subparser common arguments"""
    parser.add_argument("-o", "--outdir", default="./", help="Output diretory.")
    parser.add_argument(
        "-p",
        "--prefix",
        default="prefix",
        help="Prefix of all output files.",
    )
    parser.add_argument("-s", "--species", help="Species")
    # parser.add_argument("--thread", help="", default=4)

    return parser


def get_input(
    analysis, module, cwd, input=None, subobj=None, type=None, useColname=None
):
    input_Clustering = Path(
        os.path.join(
            cwd,
            "result",
            module,
            type,
            "Clustering",
            "seurat.h5seurat",
        )
    )
    input_Manualanno = Path(
        os.path.join(
            cwd,
            "result",
            module,
            type,
            "Manualanno",
            "seurat.h5seurat",
        )
    )
    input_QC = Path(
        os.path.join(
            cwd,
            "result",
            module,
            type,
            "QC",
            "seurat.h5seurat",
        )
    )
    if module == "M1":
        if input == None:
            if analysis == "Clustering" and subobj != None:
                if useColname == None or useColname != "clusters":
                    if input_Manualanno.exists():
                        return input_Manualanno
                    else:
                        return input_Clustering
                else:
                    return input_Clustering
            elif analysis == "Clustering" and subobj == None:
                return input_QC
            elif analysis == "Manualanno":
                return input_Clustering
            elif analysis == "Diffexp" and subobj != None:
                if useColname == "new_celltype" or useColname == None:
                    if input_Manualanno.exists():
                        return input_Manualanno
                    else:
                        return input_Clustering
            elif analysis == "Diffexp" and subobj == None:
                if input_Manualanno.exists():
                    return input_Manualanno
                else:
                    return input_Clustering
            elif analysis == "Enrich":
                sys.exit("Error: Please input the input file for Enrich analysis.")
            elif analysis == "Markerplot":
                if input_Manualanno.exists():
                    return input_Manualanno
                elif input_Clustering.exists():
                    return input_Clustering
                else:
                    sys.exit(
                        "Error: Please input the input file for Markerplot analysis."
                    )
        else:
            return input

    if module != "M1":
        if input == None:
            input_Manualanno = Path(f"{cwd}/result/M1/{type}/Manualanno/data_ob_v3.rds")
            input_Clustering = Path(f"{cwd}/result/M1/{type}/Clustering/data_ob_v3.rds")
            if input_Manualanno.exists():
                return input_Manualanno
            elif input_Clustering.exists():
                print(
                    f"Since {input_Manualanno} does not exist, {input_Clustering} is used"
                )
                return input_Clustering
            else:
                sys.exit(f"Error: please input the input file for {analysis} analysis.")
        else:
            return input


class Step:
    """
    Step class
    """

    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
        self.prefix = args.prefix
        self.species = args.species
        # self.assay = args.subparser_assay
        # self.thread = int(args.thread)

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print("Bye")


class ModuleFun:
    """
    Module class
    """

    def __init__(self, data, input, analysis, type):
        self,
        self.cwd = os.getcwd()
        self.data = data
        self.input = input
        self.analysis = analysis
        self.type = type
        self.module = self.data["module"]
        self.subnew_celltype = self.data["param"]["subnew_celltype"]
        self.subsampleid = self.data["param"]["subsampleid"]
        self.subgroup = self.data["param"]["subgroup"]
        self.subcluster = self.data["param"]["subcluster"]
        self.ncores = NCORES[self.analysis]
        self.script = SCRIPT[self.analysis]
        self.environment = ENVIRONMENT[self.analysis]
        self.category = "/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls"
        self.outdir = os.path.join(self.cwd, "result", self.data["module"])
        self.prefix = self.data["param"]["prefix"]
        self.species = self.data["species"]
        self.programID = self.data["programID"]
        try:
            self.geneanno = GENOME["geneanno"][self.data["species"]]
            self.cellanno = GENOME["cellanno"][self.data["species"]]
            self.go_bg = GENOME["go_bg"][self.data["species"]]
            self.kegg_bg = GENOME["kegg_bg"][self.data["species"]]
            self.ppispecies = GENOME["ppispecies"][self.data["species"]]
            self.gtf = GENOME["gtf"][self.data["species"]]
        except KeyError as e:
            sys.exit(
                f"Error: Please choose one of [human,mouse,rat,dog,mmul(猕猴),macaca_fascicularis(食蟹猴),sus_scrofa(猪),gallus_gallus(鸡),oar(绵羊),Capra_hircus(山羊),Caenorhabditis_elegans(秀丽隐杆线虫),silk(家蚕),Brara_Chiifu(白菜)], and the author(chenglong.liu@oebiotech.com) can be contacted to add {e}."
            )

        if self.module != "M1":
            self.outdir = os.path.join(self.outdir, self.analysis, self.type)

            if self.subnew_celltype == None:
                self.subnew_celltype = "all"
            else:
                self.outdir = os.path.join(
                    self.outdir, self.subnew_celltype.replace(",", "_")
                )
                self.type = self.type + "_" + self.subnew_celltype.replace(",", "_")
            if self.subsampleid == None:
                self.subsampleid = "all"
            else:
                self.outdir = os.path.join(
                    self.outdir, self.subsampleid.replace(",", "_")
                )
                self.type = self.type + "_" + self.subsampleid.replace(",", "_")
            if self.subgroup == None:
                self.subgroup = "all"
            else:
                self.outdir = os.path.join(self.outdir, self.subgroup.replace(",", "_"))
                self.type = self.type + "_" + self.subgroup.replace(",", "_")

            if self.subcluster == None:
                self.subcluster = "all"
            else:
                self.outdir = os.path.join(
                    self.outdir, self.subcluster.replace(",", "_")
                )
                self.type = self.type + "_" + self.subcluster.replace(",", "_")

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    def save_script(self, shell_script_content):
        shell_path = f"{self.cwd}/script/{self.module}/{self.analysis}_{self.prefix}.sh"
        with open(
            shell_path,
            "w",
        ) as file:
            file.write(shell_script_content)

        string = f"""echo qsub -V -cwd -pe  smp {self.ncores} -q scrna -o {self.cwd}/script/log/ -e {self.cwd}/script/log/ {shell_path} >> {self.cwd}/script/qsub.sh"""
        os.system(string)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print(f"Init shell script for {self.analysis} of {self.type} finished.")
