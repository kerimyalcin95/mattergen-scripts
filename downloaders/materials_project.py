from .base import BaseDownloader
from .utils import save_structure, create_metadata_row

from pathlib import Path
import pandas as pd
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
    ) -> pd.DataFrame:
        """
        Download all structures for a chemical system from the
        Materials Project.
        """

        if not self.api_key:
            raise RuntimeError(
                "Materials Project requires --api-key."
            )

        metadata = []

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

                    metadata.append(
                        create_metadata_row(
                            database="MaterialsProject",
                            chemsys=chemsys,
                            source_id=str(document.material_id),
                            filename=filename,
                            structure=document.structure,
                            formula=document.formula_pretty,
                            elements=document.elements,
                            space_group=document.symmetry.symbol,
                            space_group_number=document.symmetry.number,
                            crystal_system=document.symmetry.crystal_system,
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

        return metadata_df
