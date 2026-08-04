import pandas as pd
from pathlib import Path
from abc import ABC, abstractmethod

class BaseDownloader(ABC):

    @abstractmethod
    def download(
        self,
        chemsys: str,
        output_folder: Path,
    ) -> pd.DataFrame:
        pass