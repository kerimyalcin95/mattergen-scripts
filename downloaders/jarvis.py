from .optimade import OPTIMADEDownloader


class JARVISDownloader(OPTIMADEDownloader):

    database_name = "JARVIS"
    provider_name = "jarvis"

    BASE_URL = "https://jarvis.nist.gov/optimade/jarvisdft/v1"

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:

        energy_above_hull = attributes.get("_jarvis_ehull")

        return {
            "band_gap": attributes.get("_jarvis_band_gap"),
            "formation_energy_per_atom": attributes.get(
                "_jarvis_formation_energy_peratom"
            ),
            "energy_per_atom": attributes.get(
                "_jarvis_energy_peratom"
            ),
            "energy_above_hull": energy_above_hull,
            "is_stable": (
                energy_above_hull == 0
                if energy_above_hull is not None
                else None
            ),
        }
