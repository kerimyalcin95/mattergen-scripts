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

from downloaders.materials_project import MaterialsProjectDownloader

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

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


DOWNLOADERS = {
    "mp": MaterialsProjectDownloader
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

        if args.database == "optimade":
            kwargs["provider"] = args.provider

        downloader = DOWNLOADERS[args.database](
            api_key=args.api_key,
        )
        downloader.download(
            chemsys=chemsys,
            output_folder=output_folder,
        )

    print(f"Finished {chemsys}")

    print()

    print("All downloads finished.")


if __name__ == "__main__":
    main()
