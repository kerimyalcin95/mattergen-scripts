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

from downloaders.materials_project import MaterialsProjectDownloader
from downloaders.utils import write_metadata

import argparse
from pathlib import Path

import os

from dotenv import load_dotenv

load_dotenv()

DOWNLOADERS = {
    "mp": MaterialsProjectDownloader,
}


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Download structures from crystal databases."
    )

    parser.add_argument(
        "--database",
        required=True,
        choices=DOWNLOADERS.keys(),
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


def main():

    args = parse_arguments()

    downloader = DOWNLOADERS[args.database](
        api_key=args.api_key,
    )

    for chemsys in args.chemsys:

        print()

        print("=" * 70)
        print(f"Database : {downloader.database_name}")
        print(f"Chemsys  : {chemsys}")
        print("=" * 70)

        output_folder = downloader.create_output_folder(
            chemsys=chemsys,
            output_root=args.output,
        )

        metadata = downloader.download(
            chemsys=chemsys,
            output_folder=output_folder,
        )

        write_metadata(
            metadata,
            output_folder,
        )

    print(f"Finished {chemsys}")

    print()

    print("All downloads finished.")


if __name__ == "__main__":
    main()
