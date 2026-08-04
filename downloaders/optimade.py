from __future__ import annotations

from abc import abstractmethod
from pathlib import Path

import pandas as pd
import requests
from pymatgen.core import Lattice, Structure

from .base import BaseDownloader
from .utils import (
    create_metadata_row,
    save_structure,
)


class OPTIMADEDownloader(BaseDownloader):
    """
    Generic downloader for OPTIMADE providers.
    """

    BASE_URL: str = ""

    PAGE_SIZE = 100

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if not self.BASE_URL:
            raise ValueError(
                "BASE_URL must be defined by subclasses."
            )

    def build_filter(
        self,
        chemsys: str,
    ) -> str:
        """
        Build an OPTIMADE filter from a chemical system.

        Example
        -------
        Al-O -> elements HAS ALL "Al", "O"
        """

        elements = chemsys.split("-")

        quoted = ", ".join(
            f'"{element}"'
            for element in elements
        )

        return f"elements HAS ALL [{quoted}]"

    def iter_entries(
        self,
        chemsys: str,
    ):
        """
        Iterate over all matching OPTIMADE structures.
        """

        url = f"{self.BASE_URL}/structures"

        params = {
            "filter": self.build_filter(chemsys),
            "page_limit": self.PAGE_SIZE,
        }

        while url:

            response = requests.get(
                url,
                params=params,
                timeout=120,
            )

            response.raise_for_status()

            payload = response.json()

            for entry in payload.get("data", []):
                yield entry

            next_link = payload.get("links", {}).get("next")

            if isinstance(next_link, dict):
                url = next_link.get("href")
            else:
                url = next_link

            params = None

    def build_structure(
        self,
        attributes: dict,
    ) -> Structure:
        """
        Convert an OPTIMADE structure into a pymatgen Structure.
        """

        lattice = Lattice(
            attributes["lattice_vectors"]
        )

        return Structure(
            lattice=lattice,
            species=attributes["species_at_sites"],
            coords=attributes["cartesian_site_positions"],
            coords_are_cartesian=True,
        )

    def download(
        self,
        chemsys: str,
        output_folder: Path,
    ) -> pd.DataFrame:

        metadata = []

        downloaded = 0
        skipped = 0
        failed = 0

        for entry in self.iter_entries(chemsys):

            try:

                attributes = entry["attributes"]

                structure = self.build_structure(
                    attributes
                )

                filename = f"{entry['id']}.cif"

                filepath = output_folder / filename

                if filepath.exists():
                    skipped += 1
                else:
                    save_structure(
                        structure,
                        filepath,
                    )
                    downloaded += 1

                provider_metadata = self.extract_provider_metadata(
                    attributes
                )

                metadata.append(
                    create_metadata_row(
                        database=self.database_name,
                        chemsys=chemsys,
                        source_id=entry["id"],
                        filename=filename,
                        structure=structure,
                        formula=attributes.get(
                            "chemical_formula_reduced"
                        ),
                        elements=attributes.get(
                            "elements",
                            [],
                        ),
                        **provider_metadata,
                    )
                )

            except Exception as exception:

                failed += 1

                print(
                    f"Failed {entry.get('id')}: "
                    f"{exception}"
                )

        metadata_df = pd.DataFrame(metadata)

        if not metadata_df.empty:

            metadata_df.sort_values(
                by="source_id",
                inplace=True,
            )

        print()
        print("=" * 70)
        print(f"{self.database_name} Download Summary")
        print("=" * 70)
        print(f"Chemical system : {chemsys}")
        print(f"Downloaded      : {downloaded}")
        print(f"Skipped         : {skipped}")
        print(f"Failed          : {failed}")
        print(f"Metadata rows   : {len(metadata_df)}")
        print("=" * 70)

        return metadata_df

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:
        """
        Provider-specific metadata.

        Override in subclasses.
        """

        return {}
