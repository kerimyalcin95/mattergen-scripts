import pandas as pd
from pathlib import Path
from abc import ABC, abstractmethod


class BaseDownloader(ABC):

    database_name: str

    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def download(
        self,
        chemsys: str,
        output_folder: Path,
    ) -> pd.DataFrame:
        pass

    def create_output_folder(
        self,
        chemsys: str,
        output_root: Path,
    ) -> Path:

        folder = output_root / self.database_name / chemsys

        folder.mkdir(
            parents=True,
            exist_ok=True,
        )

        return folder
