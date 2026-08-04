from .optimade import OPTIMADEDownloader


class OQMDDownloader(OPTIMADEDownloader):

    database_name = "OQMD"
    provider_name = "oqmd"

    BASE_URL = "https://oqmd.org/optimade/v1"

    def extract_provider_metadata(
        self,
        attributes: dict,
    ) -> dict:

        return {
            "band_gap": attributes.get("_oqmd_band_gap"),
            "energy_above_hull": attributes.get("_oqmd_stability"),
            "is_stable": attributes.get("_oqmd_stable"),
        }
