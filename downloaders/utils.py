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
