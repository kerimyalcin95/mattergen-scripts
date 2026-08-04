from .optimade import OPTIMADEDownloader


class AlexandriaDownloader(OPTIMADEDownloader):

    database_name = "Alexandria"

    BASE_URL = (
        "https://alexandria.icams.rub.de/pbe"
    )

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:

        energy_above_hull = attributes.get(
            "_alexandria_energy_above_hull"
        )

        return {
            "band_gap": attributes.get(
                "_alexandria_band_gap"
            ),
            "formation_energy_per_atom": attributes.get(
                "_alexandria_formation_energy_per_atom"
            ),
            "energy_per_atom": attributes.get(
                "_alexandria_energy_per_atom"
            ),
            "energy_above_hull": energy_above_hull,
            "is_stable": (
                attributes.get("_alexandria_is_stable")
                if attributes.get("_alexandria_is_stable") is not None
                else (
                    energy_above_hull == 0
                    if energy_above_hull is not None
                    else None
                )
            ),
        }
