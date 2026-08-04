from __future__ import annotations

from pathlib import Path
import requests

from pymatgen.core import (
    Lattice,
    Structure,
)

from tqdm import tqdm

from .base import BaseDownloader
from .utils import (
    save_structure,
)


class OPTIMADEDownloader(BaseDownloader):
    """
    Generic downloader for OPTIMADE providers.
    """

    database_name = ""
    provider_name = ""

    BASE_URL = ""

    PAGE_SIZE = 100

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self.session = requests.Session()

        self.session.headers.update(
            {
                "User-Agent": (
                    "Mozilla/5.0 "
                    "(Windows NT 10.0; Win64; x64)"
                ),
                "Accept": "application/json",
            }
        )

        if not self.BASE_URL:
            raise ValueError(
                "BASE_URL must be defined."
            )

    ####################################################################
    # Query
    ####################################################################

    def build_filter(
        self,
        chemsys: str,
    ) -> str:

        elements = chemsys.split("-")

        quoted = ", ".join(
            f'"{element}"'
            for element in elements
        )

        return (
            f"elements HAS ALL {quoted} "
            f"AND nelements<={len(elements)}"
        )

    def iter_entries(
        self,
        chemsys: str,
    ):

        url = f"{self.BASE_URL}/structures"

        params = {
            "filter": self.build_filter(
                chemsys
            ),
            "page_limit": self.PAGE_SIZE,
        }

        while url:

            response = self.session.get(
                url,
                params=params,
                timeout=120,
            )

            print("STATUS :", response.status_code)
            print("FINAL  :", response.url)

            if response.status_code != 200:
                print(
                    f""
                    f"{self.database_name} returned "
                    f"{response.status_code}: {response.url}"
                )
                return

            payload = response.json()

            for entry in payload.get(
                "data",
                [],
            ):
                yield entry

            next_link = (
                payload
                .get("links", {})
                .get("next")
            )

            if isinstance(next_link, dict):
                url = next_link.get("href")
            else:
                url = next_link

            params = None

    ####################################################################
    # Structure
    ####################################################################

    def has_complete_structure(
        self,
        attributes: dict,
    ) -> bool:

        return (
            attributes.get(
                "lattice_vectors"
            )
            and (
                attributes.get(
                    "cartesian_site_positions"
                )
                or attributes.get(
                    "fractional_site_positions"
                )
            )
            and attributes.get(
                "species_at_sites"
            )
        )

    def build_structure(
        self,
        attributes: dict,
    ) -> Structure:

        lattice = Lattice(
            attributes["lattice_vectors"]
        )

        coords = attributes.get(
            "cartesian_site_positions"
        )

        coords_are_cartesian = True

        if coords is None:

            coords = attributes[
                "fractional_site_positions"
            ]

            coords_are_cartesian = False

        return Structure(
            lattice=lattice,
            species=attributes[
                "species_at_sites"
            ],
            coords=coords,
            coords_are_cartesian=coords_are_cartesian,
        )

    def download_structure(
        self,
        entry: dict,
    ) -> Structure | None:

        attributes = entry["attributes"]

        if not self.has_complete_structure(
            attributes
        ):
            return None

        return self.build_structure(
            attributes
        )

    ####################################################################
    # Metadata
    ####################################################################

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:

        return {}

    ####################################################################
    # Download
    ####################################################################

    def download(
        self,
        chemsys: str,
        output_folder: Path,
    ) -> None:

        downloaded = 0
        skipped = 0
        failed = 0

        print(
            f"Connecting to {self.database_name}..."
        )

        for entry in tqdm(
            self.iter_entries(
                chemsys,
            ),
            desc=chemsys,
        ):

            try:

                attributes = entry[
                    "attributes"
                ]

                structure = (
                    self.download_structure(
                        entry
                    )
                )

                if structure is None:
                    continue

                filename = (
                    f"{entry['id']}.cif"
                )

                filepath = (
                    output_folder
                    / filename
                )

                if filepath.exists():

                    skipped += 1

                else:

                    save_structure(
                        structure,
                        filepath,
                    )

                    downloaded += 1

            except Exception as exception:

                failed += 1

                print(
                    f"Failed "
                    f"{entry.get('id')}: "
                    f"{exception}"
                )

        print()
        print("=" * 70)
        print(
            f"{self.database_name} Download Summary"
        )
        print("=" * 70)
        print(
            f"Chemical system : {chemsys}"
        )
        print(
            f"Downloaded      : {downloaded}"
        )
        print(
            f"Skipped         : {skipped}"
        )
        print(
            f"Failed          : {failed}"
        )
        print("=" * 70)