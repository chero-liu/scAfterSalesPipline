#!/usr/bin/env python3
# encoding: utf-8
"""
Author  : chenglong liu
Contact : chenglong.liu@oebiotech.com
File   : ppi.py
Time   : 202405
"""

from typing import Optional
import pandas as pd
from tools.utils import check_input


@check_input
def load_diffgene_data(
    input: Optional[str],
    noft: Optional[int] = None,
    gene_col: Optional[str] = None,
    fold_change_col: Optional[str] = None,
    regulation_col: Optional[str] = None,
    up_regulation_value: Optional[str] = None,
    down_regulation_value: Optional[str] = None,
) -> pd.DataFrame:
    """
    Load differential gene expression data and select top records by fold change.

    Parameters:
    - input: str
        Path to the file containing differential gene expression data.
    - noft: Optional[int] (default None)
        Number of top records to select based on the fold change column for each regulation type.
        If None, all records are returned without filtering.
    - gene_col: str (default 'gene')
        The column name for genes in the dataset.
    - fold_change_col: str (default 'FoldChange')
        The column name for fold change values in the dataset.
    - regulation_col: str (default 'Regulation')
        The column name for regulation types in the dataset.
    - up_regulation_value: str (default 'Up')
        The value that represents up-regulation in the dataset.
    - down_regulation_value: str (default 'Down')
        The value that represents down-regulation in the dataset.

    Returns:
    - pd.DataFrame
        DataFrame containing the differential gene expression data, possibly filtered by the top fold changes.
    """
    # Load the data
    diffgene = pd.read_table(input)

    # Keep only relevant columns
    try:
        diffgene = diffgene[[gene_col, fold_change_col, regulation_col]]
    except KeyError as e:
        raise KeyError(f"One or more specified columns do not exist in the data: {e}")

    # If noft is specified, filter the top fold changes
    if noft is not None:
        # Check that noft is a positive integer
        if not isinstance(noft, int) or noft <= 0:
            raise ValueError("noft must be a positive integer")

        # Split and filter the dataframe by regulation type
        up_reg = diffgene[diffgene[regulation_col] == up_regulation_value]
        down_reg = diffgene[diffgene[regulation_col] == down_regulation_value]

        try:
            up = up_reg.nlargest(noft, fold_change_col)
            down = down_reg.nsmallest(noft, fold_change_col)
        except KeyError as e:
            raise KeyError(f"Column for sorting does not exist in the data: {e}")

        diffgene = pd.concat([up, down], ignore_index=True)

    return diffgene


@check_input
def load_stringdb(
    input: Optional[str],
    sep: str = "\t",
    header: Optional[int] = 1,
    protein1_col: Optional[str] = "protein1",
    protein2_col: Optional[str] = "protein2",
    combined_score_col: Optional[str] = "combined_score",
) -> pd.DataFrame:
    """Load STRING database gene interaction data for a specific species.

    Parameters:
        base_path (str): The base directory where the species folders are located.
        species (Optional[str]): The species folder name. If None, no species subfolder will be used.
        filename (str): The name of the file inside the species directory.
        sep (str): The delimiter of the data file.
        header (int): Row to use as header.
        protein1_col (str): Column name for the first protein.
        protein2_col (str): Column name for the second protein.
        combined_score_col (str): Column name for the combined score.

    Returns:
        pd.DataFrame: DataFrame containing the loaded STRING database data.
    """
    # Read the data file
    stringdb = pd.read_table(input, sep=sep, header=header)

    # Fix the header and first row if needed
    if header is None:
        colnames = [protein1_col, protein2_col, combined_score_col]
        stringdb.columns = colnames
    else:
        # Rename the default columns to given parameters
        rename_columns = {
            stringdb.columns[0]: protein1_col,
            stringdb.columns[1]: protein2_col,
            stringdb.columns[2]: combined_score_col,
        }
        stringdb.rename(columns=rename_columns, inplace=True)

    return stringdb
