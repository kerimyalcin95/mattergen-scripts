from .base import BaseDownloader
from .utils import save_structure

from pathlib import Path
from tqdm import tqdm

# Materials Project
from mp_api.client import MPRester

import os

from dotenv import load_dotenv

load_dotenv()


class MaterialsProjectDownloader(BaseDownloader):

    database_name = "MaterialsProject"

    def __init__(self, api_key: str):
        self.api_key = api_key

    def download(
        self,
        chemsys: str,
        output_folder: Path,
    ) -> None:
        """
        Download all structures for a chemical system from the
        Materials Project.
        """

        if not self.api_key:
            raise RuntimeError(
                "Materials Project requires --api-key."
            )

        downloaded = 0
        skipped = 0
        failed = 0

        print("Connecting to Materials Project...")

        with MPRester(self.api_key) as mpr:

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

                except Exception as exception:

                    failed += 1

                    print(
                        f"Failed {document.material_id}: "
                        f"{exception}"
                    )

        print()
        print("=" * 70)
        print("Materials Project Download Summary")
        print("=" * 70)
        print(f"Chemical system : {chemsys}")
        print(f"Downloaded      : {downloaded}")
        print(f"Skipped         : {skipped}")
        print(f"Failed          : {failed}")
        print("=" * 70)
