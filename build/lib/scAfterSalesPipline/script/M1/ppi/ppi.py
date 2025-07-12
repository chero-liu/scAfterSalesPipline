#!/usr/bin/env python3
# encoding: utf-8
"""
Author  : chenglong liu
Contact : chenglong.liu@oebiotech.com
File   : ppi.py
Time   : 202405
"""

import argparse
import os
from tools.load import load_diffgene_data, load_stringdb
from tools.utils import get_network
import sys


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate PPI network from STRING database for differential gene expression data."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input xls file with differential gene."
    )
    parser.add_argument(
        "-o", "--outdir", required=True, help="Output directory for results."
    )
    parser.add_argument(
        "-p", "--prefix", required=True, help="Prefix for the output files."
    )
    parser.add_argument(
        "-s",
        "--species",
        required=True,
        help="NCBI taxon identifiers, e.g., Human (9606), Mouse (10090), Rat (10116).Or dir path of gene2gene_network.xls and TF_anno.xls.",
    )
    parser.add_argument(
        "--score_top",
        type=int,
        default=300,
        help=" ",
    )
    parser.add_argument(
        "--gene_col",
        type=str,
        default="gene",
        help="Diffgene xls gene_col name.",
    )
    parser.add_argument(
        "--fold_change_col",
        type=str,
        default="FoldChange",
        help="Diffgene xls fold_change_col name.",
    )
    parser.add_argument(
        "--regulation_col",
        type=str,
        default="Regulation",
        help="Diffgene xls regulation_col name.",
    )
    parser.add_argument(
        "--up_regulation_value",
        type=str,
        default="Up",
        help="Diffgene xls up_regulation_value value.",
    )
    parser.add_argument(
        "--down_regulation_value",
        type=str,
        default="Down",
        help="Diffgene xls down_regulation_value value.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure output directory exists
    outdir = os.path.realpath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # Process differential gene expression data
    diffgene = load_diffgene_data(
        args.input,
        noft=None,
        gene_col=args.gene_col,
        fold_change_col=args.fold_change_col,
        regulation_col=args.regulation_col,
        up_regulation_value=args.up_regulation_value,
        down_regulation_value=args.down_regulation_value,
    )
    # Load the STRING PPI network data for the specified species
    if os.path.isdir(args.species):
        stringd_filepath = f"{args.species}/gene2gene_network.xls"
        TF_filepath = f"{args.species}/TF_anno.xls"
    else:
        stringd_filepath = (
            f"/gpfs/oe-database/reference/RNA/{args.species}/gene2gene_network.xls"
        )
        TF_filepath = f"/gpfs/oe-database/reference/RNA/{args.species}/TF_anno.xls"
    stringdb = load_stringdb(stringd_filepath)

    # Generate network dataframe based on differential gene and STRING data
    network_df_top = get_network(
        diffgene,
        stringdb,
        args.score_top,
        args.gene_col,
        args.regulation_col,
    )

    network_df = get_network(
        diffgene,
        stringdb,
        None,
        args.gene_col,
        args.regulation_col,
    )

    if network_df.empty | network_df_top.empty:
        sys.exit(
            "No interaction pairs between diffgene are present in stringdb. Please check for potential reasons."
        )

    # Save the resulting network DataFrame
    output_file = os.path.join(outdir, f"{args.prefix}.ppi_network.tsv")
    network_df.to_csv(output_file, sep="\t", index=False)
    output_file_top = os.path.join(
        outdir, f"{args.prefix}.ppi_network_top{args.score_top}.tsv"
    )
    network_df_top.to_csv(output_file_top, sep="\t", index=False)
    # Generate 3D network visualization
    if os.path.exists(TF_filepath):
        os.system(
            f"bash {os.path.dirname(os.path.realpath(__file__))}/script/network_3d.sh  {outdir}  {args.prefix}.ppi_network_top{args.score_top}.tsv  {TF_filepath}"
        )
    else:
        print(f"{TF_filepath} does not exist, so the network_3d is not drawn")


if __name__ == "__main__":
    main()
