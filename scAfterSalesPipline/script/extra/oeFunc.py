import pandas as pd


def get_auto_anno_info(
    input_str: str,
    attribute: str,
    idx_data_path="/gpfs/oe-scrna/xiaojiacheng/public/celltyping_ref_updated/reference_celltype_source.csv",
) -> tuple:
    # Validate input
    if not isinstance(input_str, str) or "-" not in input_str:
        input_str = f"{input_str}-"

    try:
        name, tissue = input_str.split("-")
    except ValueError:
        raise ValueError("Input should be a string with '-', e.g. '人-肝癌'.")

    tissue = tissue.replace("组织", "")

    # Load index data
    idx_data = pd.read_csv(idx_data_path)

    # Validate attribute
    if attribute not in idx_data.columns:
        raise ValueError(
            f"Attribute {attribute} not found in the index data. Choose from {idx_data.columns[2:]}."
        )

    # Validate species
    if name not in idx_data["names"].values:
        raise ValueError(f"Species {name} not found in the index data.")

    # Filter data for the given species
    species_data = idx_data[idx_data["names"] == name]

    # Define default tissue for each species
    default_tissue = "default"

    # Get result based on tissue availability
    if tissue not in species_data["tissue"].values:
        print(
            f"Tissue {tissue} not found for species {name}, choose to {default_tissue}."
        )
        tissue = default_tissue

    result = species_data.loc[species_data["tissue"] == tissue, attribute].values[0]
    species = species_data["species"].values[0]

    return species, result


def main():
    pass


if __name__ == "__main__":
    main()
