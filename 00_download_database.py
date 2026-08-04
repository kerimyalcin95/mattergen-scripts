#!/usr/bin/env python3
"""
Download crystal structures from supported databases.

Supported
---------
- Materials Project
- COD (placeholder)
- OPTIMADE (placeholder)

Example
-------
python 00_download_database.py \
    --database mp \
    --chemsys Al-O \
    --api-key YOUR_API_KEY \
    --output ../../data
"""

from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime

import pandas as pd
from tqdm import tqdm

from pymatgen.io.cif import CifWriter

# Materials Project
from mp_api.client import MPRester

DATABASES = {
    "mp": "MaterialsProject",
    "cod": "COD",
    "optimade": "OPTIMADE",
}


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Download structures from crystal databases."
    )

    parser.add_argument(
        "--database",
        required=True,
        choices=DATABASES.keys(),
        help="Database to download from.",
    )

    parser.add_argument(
        "--chemsys",
        required=True,
        nargs="+",
        help="Chemical systems (e.g. Al-O Fe-O Y-Al-O)"
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output root folder."
    )

    parser.add_argument(
        "--api-key",
        default=None,
        help="Materials Project API key."
    )

    parser.add_argument(
        "--provider",
        default=None,
        help="OPTIMADE provider."
    )

    return parser.parse_args()


def create_output_folder(
    database: str,
    chemsys: str,
    output_root: Path,
) -> Path:

    database_name = DATABASES[database]

    if database == "optimade":
        raise NotImplementedError(
            "OPTIMADE providers will be implemented in Part 2."
        )

    root = output_root / database_name / chemsys

    root.mkdir(
        parents=True,
        exist_ok=True,
    )

    return root


def write_metadata(
    dataframe: pd.DataFrame,
    output_folder: Path,
):

    dataframe.to_csv(
        output_folder / "metadata.csv",
        index=False,
    )


def download_materials_project(
    api_key: str,
    chemsys: str,
    output_folder: Path,
):
    """
    Download structures from Materials Project.
    """

    if api_key is None:
        raise RuntimeError(
            "Materials Project requires --api-key"
        )

    metadata = []

    with MPRester(api_key) as mpr:

        docs = mpr.materials.summary.search(
            chemsys=chemsys,
            fields=[
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "density",
                "volume",
                "elements",
            ],
        )

        docs = list(docs)

        print(
            f"Found {len(docs)} structures."
        )

        for doc in tqdm(docs):

            structure = doc.structure

            filename = (
                f"{doc.material_id}.cif"
            )

            path = output_folder / filename

            CifWriter(
                structure
            ).write_file(path)

            metadata.append(
                {
                    "database": "MaterialsProject",
                    "source_id": str(doc.material_id),
                    "filename": filename,
                    "formula": doc.formula_pretty,
                    "elements": ",".join(map(str, doc.elements)),
                    "space_group": doc.symmetry.symbol,
                    "crystal_system": doc.symmetry.crystal_system,
                    "density": doc.density,
                    "volume": doc.volume,
                }
            )

    write_metadata(
        pd.DataFrame(metadata),
        output_folder,
    )


def download_cod(
    chemsys: str,
    output_folder: Path,
):
    """
    Placeholder.

    Implemented in Part 2.
    """

    raise NotImplementedError(
        "COD support will be added in Part 2."
    )


def download_optimade(
    provider: str,
    chemsys: str,
    output_folder: Path,
):
    """
    Placeholder.

    Implemented in Part 2.
    """

    raise NotImplementedError(
        "OPTIMADE support will be added in Part 2."
    )


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


DOWNLOADERS = {
    "mp": download_materials_project,
    "cod": download_cod,
    "optimade": download_optimade,
}


def main():

    args = parse_arguments()

    for chemsys in args.chemsys:

        print()

        print("=" * 70)
        print(f"Database : {args.database}")
        print(f"Chemsys  : {chemsys}")
        print("=" * 70)

        output_folder = create_output_folder(
            database=args.database,
            chemsys=chemsys,
            output_root=args.output,
        )

        if args.database == "mp":

            download_materials_project(
                api_key=args.api_key,
                chemsys=chemsys,
                output_folder=output_folder,
            )

        elif args.database == "cod":

            download_cod(
                chemsys=chemsys,
                output_folder=output_folder,
            )

        elif args.database == "optimade":

            download_optimade(
                provider=args.provider,
                chemsys=chemsys,
                output_folder=output_folder,
            )

        print(f"Finished {chemsys}")

    print()

    print("All downloads finished.")


if __name__ == "__main__":
    main()
