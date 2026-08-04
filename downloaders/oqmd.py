from .optimade import OPTIMADEDownloader


class OQMDDownloader(OPTIMADEDownloader):

    database_name = "OQMD"

    BASE_URL = "https://oqmd.org/optimade/v1"