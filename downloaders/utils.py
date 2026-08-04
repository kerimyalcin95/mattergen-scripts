from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pathlib import Path
import pandas as pd


def save_structure(
    structure,
    path: Path,
):
    """
    Save a pymatgen Structure as CIF.
    """

    writer = CifWriter(
        structure,
        symprec=0.1,
    )

    writer.write_file(path)


def write_metadata(
    dataframe: pd.DataFrame,
    output_folder: Path,
):

    dataframe.to_csv(
        output_folder / "metadata.csv",
        index=False,
    )


def create_metadata_row(
    database: str,
    chemsys: str,
    source_id: str,
    filename: str,
    structure,
    formula: str,
    elements,
    symmetry,
    density,
    volume,
    **provider_metadata,
) -> dict:
    """
    Create a standardized metadata row for every database.
    """

    row = {
        "database": database,
        "chemsys": chemsys,
        "source_id": source_id,
        "filename": filename,
        "formula": formula,
        "elements": ",".join(map(str, elements)),
        "num_sites": len(structure),
        "space_group": symmetry.symbol,
        "space_group_number": symmetry.number,
        "crystal_system": symmetry.crystal_system,
        "density": density,
        "volume": volume,
    }

    row.update(provider_metadata)

    return row
