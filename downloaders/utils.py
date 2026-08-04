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
    density=None,
    volume=None,
    band_gap=None,
    energy_above_hull=None,
    is_stable=None,
    space_group=None,
    space_group_number=None,
    crystal_system=None,
) -> dict:
    """
    Create a standardized metadata row.

    If symmetry information is not provided, it is computed from the
    reconstructed structure.
    """

    if (
        space_group is None
        or space_group_number is None
        or crystal_system is None
    ):

        try:

            analyzer = SpacegroupAnalyzer(structure)

            space_group = analyzer.get_space_group_symbol()
            space_group_number = analyzer.get_space_group_number()
            crystal_system = analyzer.get_crystal_system()

        except Exception:

            space_group = None
            space_group_number = None
            crystal_system = None

    if density is None:
        density = float(structure.density)

    if volume is None:
        volume = structure.volume

    return {
        "database": database,
        "chemsys": chemsys,
        "source_id": source_id,
        "filename": filename,
        "formula": formula,
        "elements": ",".join(map(str, elements)),
        "num_sites": len(structure),
        "space_group": space_group,
        "space_group_number": space_group_number,
        "crystal_system": crystal_system,
        "density": density,
        "volume": volume,
        "band_gap": band_gap,
        "energy_above_hull": energy_above_hull,
        "is_stable": is_stable,
    }
