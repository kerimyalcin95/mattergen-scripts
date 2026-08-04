from pathlib import Path
import tempfile

from pymatgen.core import Structure

from .optimade import OPTIMADEDownloader


class CODDownloader(OPTIMADEDownloader):

    database_name = "COD"
    provider_name = "cod"

    BASE_URL = "https://www.crystallography.net/cod/optimade"

    CIF_URL = (
        "https://www.crystallography.net/cod/{cod_id}.cif"
    )

    def download_structure(
        self,
        entry: dict,
    ) -> Structure | None:
        """
        Download the original CIF from COD and construct the
        pymatgen Structure from it.

        The COD OPTIMADE endpoint only provides metadata and is
        therefore used only for searching.
        """

        attributes = entry["attributes"]

        cod_id = attributes.get(
            "_cod_file",
            entry["id"],
        )

        response = self.session.get(
            self.CIF_URL.format(
                cod_id=cod_id,
            ),
            timeout=120,
        )

        response.raise_for_status()

        with tempfile.NamedTemporaryFile(
            suffix=".cif",
            delete=False,
        ) as temporary:

            temporary.write(response.content)

            temporary_path = Path(
                temporary.name
            )

        try:

            return Structure.from_file(
                temporary_path
            )

        finally:

            temporary_path.unlink(
                missing_ok=True,
            )

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:

        return {}