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

import pandas as pd
from tqdm import tqdm

from pymatgen.io.cif import CifWriter

# Materials Project
from mp_api.client import MPRester

import os

from dotenv import load_dotenv

load_dotenv()

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
        default=os.getenv("MP_API_KEY"),
        help="Defaults to MP_API_KEY from .env",
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
    Download all structures for a chemical system from the
    Materials Project.
    """

    if not api_key:
        raise RuntimeError(
            "Materials Project requires --api-key."
        )

    metadata = []

    downloaded = 0
    skipped = 0
    failed = 0

    print("Connecting to Materials Project...")

    with MPRester(api_key) as mpr:

        documents = mpr.materials.summary.search(
            chemsys=chemsys,
            fields=[
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "density",
                "volume",
                "elements",
                "band_gap",
                "energy_above_hull",
                "is_stable",
            ],
        )

        for document in tqdm(
            documents,
            desc=f"{chemsys}",
        ):

            filename = f"{document.material_id}.cif"
            filepath = output_folder / filename

            try:

                if filepath.exists():
                    skipped += 1
                else:

                    save_structure(
                        document.structure,
                        filepath,
                    )

                    downloaded += 1

                symmetry = document.symmetry

                metadata.append(
                    create_metadata_row(
                        database="MaterialsProject",
                        chemsys=chemsys,
                        source_id=str(document.material_id),
                        filename=filename,
                        structure=document.structure,
                        formula=document.formula_pretty,
                        elements=document.elements,
                        symmetry=symmetry,
                        density=document.density,
                        volume=document.volume,
                        band_gap=getattr(document, "band_gap", None),
                        energy_above_hull=getattr(
                            document, "energy_above_hull", None),
                        is_stable=getattr(document, "is_stable", None),
                    )
                )

            except Exception as exception:

                failed += 1

                print(
                    f"Failed {document.material_id}: "
                    f"{exception}"
                )

    metadata_df = pd.DataFrame(metadata)

    metadata_df.sort_values(
        by="source_id",
        inplace=True,
    )

    write_metadata(
        metadata_df,
        output_folder,
    )

    print()
    print("=" * 70)
    print("Materials Project Download Summary")
    print("=" * 70)
    print(f"Chemical system : {chemsys}")
    print(f"Downloaded      : {downloaded}")
    print(f"Skipped         : {skipped}")
    print(f"Failed          : {failed}")
    print(f"Metadata rows   : {len(metadata_df)}")
    print("=" * 70)


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
    band_gap=None,
    energy_above_hull=None,
    is_stable=None,
) -> dict:
    """
    Create a standardized metadata row for every database.
    """

    return {
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
        "band_gap": band_gap,
        "energy_above_hull": energy_above_hull,
        "is_stable": is_stable,
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

        kwargs = {
            "chemsys": chemsys,
            "output_folder": output_folder,
        }

        if args.database == "mp":
            kwargs["api_key"] = args.api_key

        if args.database == "optimade":
            kwargs["provider"] = args.provider

        DOWNLOADERS[args.database](**kwargs)

        print(f"Finished {chemsys}")

    print()

    print("All downloads finished.")


if __name__ == "__main__":
    main()
