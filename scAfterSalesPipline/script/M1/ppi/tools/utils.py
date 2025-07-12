#!/usr/bin/env python3
# encoding: utf-8
"""
Author  : chenglong liu
Contact : chenglong.liu@oebiotech.com
File   : ppi.py
Time   : 202405
"""

import os
import sys
import pandas as pd
from typing import Optional, Callable
from functools import wraps


def get_network(
    diffgene: pd.DataFrame,
    stringdb: pd.DataFrame,
    top: Optional[int],
    gene_col: Optional[str],
    regulation_col: Optional[str],
) -> pd.DataFrame:
    """Construct a network DataFrame by cross-referencing STRING database with differential gene data."""
    # Input validation
    required_columns_diffgene = {gene_col, regulation_col}
    required_columns_stringdb = {"protein1", "protein2", "combined_score"}

    if not required_columns_diffgene.issubset(diffgene.columns):
        missing = required_columns_diffgene - set(diffgene.columns)
        raise ValueError(f"diffgene DataFrame is missing columns: {missing}")

    if not required_columns_stringdb.issubset(stringdb.columns):
        missing = required_columns_stringdb - set(stringdb.columns)
        raise ValueError(f"stringdb DataFrame is missing columns: {missing}")

    # Avoid using inplace modifications to ensure DataFrame integrity
    # Map regulations to the corresponding genes in the stringdb dataframe
    regulation_dict = pd.Series(
        diffgene[regulation_col].values, index=diffgene[gene_col]
    ).to_dict()

    # Filter interactions where both proteins are in the diffgene dataset
    network = (
        stringdb[
            stringdb["protein1"].isin(diffgene[gene_col])
            & stringdb["protein2"].isin(diffgene[gene_col])
        ]
        .assign(
            Regulation_A=lambda x: x["protein1"].map(regulation_dict),
            Regulation_B=lambda x: x["protein2"].map(regulation_dict),
        )
        .rename(
            columns={
                "protein1": "preferredName_A",
                "protein2": "preferredName_B",
                "combined_score": "score",
            }
        )[
            [
                "preferredName_A",
                "preferredName_B",
                "Regulation_A",
                "Regulation_B",
                "score",
            ]
        ]
    )

    # Ensure that score can be converted to float and normailze it
    try:
        network["score"] = network["score"].astype(float) / 1000
    except ValueError:
        raise ValueError("Score column could not be converted to float")

    # Filter top interactions if specified
    if top is not None:
        network = network.nlargest(top, "score")

    return network


def check_input(func: Callable) -> Callable:
    @wraps(func)
    def wrapper_check_input(input: Optional[str], *args, **kwargs) -> None:
        """
        Check if the given path is a directory or a file path.

        If it is a directory path and it does not exist, create the directory.
        If it is a file path, check that the file exists or not.
        If the file does not exist, issue a warning and exit the program.

        Parameters:
        - input (str): The path to check.

        Returns:
        None
        """
        if input is None:
            sys.exit("No input path provided.")

        # Get the absolute path
        absolute_path = os.path.realpath(input)

        if not os.path.exists(absolute_path):
            if os.path.splitext(absolute_path)[1]:  # Path has a file extension
                sys.exit(f"The file at path {absolute_path} does not exist.")
            else:
                # Assume it's a directory path, create the directory
                os.makedirs(absolute_path, exist_ok=True)
                print(f"Directory created at {absolute_path}")
        else:
            if os.path.isdir(absolute_path):
                print(f"Directory already exists at {absolute_path}")
            else:
                print(f"File exists at {absolute_path}")

        # Call the decorated function if checks pass
        return func(absolute_path, *args, **kwargs)

    return wrapper_check_input
