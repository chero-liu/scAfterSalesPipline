import glob
import os
import shutil
import unittest
from pathlib import Path
import pandas as pd
from scAfterSalesPipline.tools.utils import Step, s_common


class Clean(Step):
    def __init__(self, args):
        super().__init__(args)
        self.input = args.input + "/" if args.input != "./" else args.input

    def check_file_exists(self, target_file):
        for dirpath, dirnames, filenames in os.walk(self.input):
            if target_file in filenames:
                print(f"Warning: The file '{target_file}' exists at: {dirpath}")
                return

    def security_measure(self):
        """
        This function copies all files and directories from the current working directory
        to a new directory named 'copy' located at the same level as the current directory,
        provided the input is one of the specified trigger strings.
        """
        trigger_conditions = {".", "./", ".//"}

        if self.input in trigger_conditions:
            current_directory = os.getcwd()
            parent_directory = os.path.dirname(current_directory)
            target_directory = os.path.join(parent_directory, "copy")

            os.makedirs(target_directory, exist_ok=True)

            for item in os.listdir(current_directory):
                source_item = os.path.join(current_directory, item)
                destination_item = os.path.join(target_directory, item)
                if os.path.isfile(source_item):
                    shutil.copy2(source_item, destination_item)
                elif os.path.isdir(source_item):
                    shutil.copytree(source_item, destination_item, dirs_exist_ok=True)
            self.input = target_directory
            print(f"All items have been copied to: {target_directory}")

    def delete_specific_file_in_directory(self, dir_name):
        for root, dirs, files in os.walk(self.input):
            if (
                dir_name in dirs
                and "harmony_Dimension_Reduction" not in dirs
                and Path(f"{root}/sampleid-batchid.xls").exists()
            ):

                os.system(f"rm {root}/sampleid-batchid.xls")

    def run(self):
        files_to_delete = [
            # ------------Clustering----------------
            "*.h5seurat",
            "*.rds",
            "*.h5ad",
            "mnn_Dimension_Reduction",
            "pca_Dimension_Reduction",
            "harmony_Dimension_Reduction",
            "output.json.tsv",
            # ------------Findmarker----------------
            "all_markers_for_each_cluster.xls",
            "top10_markers_for_each_cluster.xls",
            # ------------Diffenrich----------------
            "enrichment_sh",
            # ------------Diffexp----------------
            "*-all_diffexp_genes.xls",
            "*-diff-pval-*-FC-*[0-9].xls",
            "*-diff-pvalue-*-FC-*[0-9].xls",
            "*-diff-padj-*-FC-*[0-9].xls",
            # ------------log----------------
            "Rplots.pdf",
            "*.sh*",
            "*.sh",
            # ------------GSVA----------------
            "GSVA_enrichment_results.xls",
            # ------------Enrich----------------
            "enrichment-*-gene_module_*-Total.circos.pdf",
            "enrichment-*-cluster_*-Total.circos.pdf",
        ]

        self.security_measure()

        # clean

        self.delete_specific_file_in_directory("pca_Dimension_Reduction")

        for pattern in files_to_delete:
            file_paths = glob.glob(
                os.path.join(self.input, "**", pattern), recursive=True
            )
            for file_path in file_paths:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    print(f"Deleted file: {file_path}")
                else:
                    shutil.rmtree(file_path)
                    print(f"Deleted directory: {file_path}")

        # os.system(f"find {self.input}* -empty -delete")
        os.system(f"find {self.input}* -empty -exec echo {{}} \; -delete")

        # check
        print("Checking")
        diffenrich_expected_count = 23
        diffenrich_results_stat_paths = glob.glob(
            self.input + "/**/*_enrichment/group_*/", recursive=True
        )
        for diffenrich_results_stat_path in diffenrich_results_stat_paths:
            files = os.listdir(diffenrich_results_stat_path)
            file_count = len(files)

            if file_count != diffenrich_expected_count:
                print(
                    f"Warning: The folder '{diffenrich_results_stat_path}' contains {file_count} files. Expected {diffenrich_expected_count} files."
                )

        diffgene_conuts = 20
        diffexp_results_stat_paths = glob.glob(
            self.input + "/**/diffexp_results_stat.xls", recursive=True
        )
        for diffexp_results_stat_path in diffexp_results_stat_paths:
            diffexp_results_stat = pd.read_csv(diffexp_results_stat_path, sep="\t")
            if diffexp_results_stat.shape[0] != 1:
                print(f"Warning: The file '{diffexp_results_stat_path}' rows != 1.")
            if diffexp_results_stat.iloc[0, 4] <= diffgene_conuts:
                print(
                    f"Warning: The file '{diffexp_results_stat_path}' diff genes <= {diffgene_conuts}."
                )

        diffppi_conuts = 9
        folders = glob.glob(self.input + "/**/*ppi*/*/", recursive=True)
        for folder in folders:
            files = os.listdir(folder)
            file_count = len(files)

            if file_count != diffppi_conuts:
                print(
                    f"Warning: The folder '{folder}' contains {file_count} files. Expected {diffppi_conuts} files."
                )

        self.check_file_exists("pdf_overlap_checker.log")


def clean(args):
    with Clean(args) as runner:
        runner.run()


def get_opts_clean(parser, sub_program=True):
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
