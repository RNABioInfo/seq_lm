import pandas as pd

from dataclasses import dataclass
from pathlib import Path

@dataclass
class FlagstatResult:
    total_reads: int
    primary_reads: int
    secondary_reads: int
    supplementary_reads: int
    total_mapped: int
    primary_mapped: int

    @property
    def total_mapped_fraction(self) -> float:
        return self.total_mapped / self.total_reads

    @property
    def primary_mapped_fraction(self) -> float:
        return self.primary_mapped / self.primary_reads

    @classmethod
    def from_df(cls, df: pd.DataFrame) -> "FlagstatResult":
        return cls(
            total_reads=int(df.iat[0, 0]),  # type: ignore
            primary_reads=int(df.iat[1, 0]),  # type: ignore
            secondary_reads=int(df.iat[2, 0]),  # type: ignore
            supplementary_reads=int(df.iat[3, 0]),  # type: ignore
            total_mapped=int(df.iat[6, 0]),  # type: ignore
            primary_mapped=int(df.iat[8, 0]),  # type: ignore
        )

    @classmethod
    def from_tsv(cls, tsv_path: Path) -> "FlagstatResult":
        return cls.from_df(pd.read_csv(tsv_path, sep="\t", header=None))

    def __add__(self, other: "FlagstatResult") -> "FlagstatResult":
        return FlagstatResult(
            total_reads=self.total_reads + other.total_reads,
            primary_reads=self.primary_reads + other.primary_reads,
            secondary_reads=self.secondary_reads + other.secondary_reads,
            supplementary_reads=self.supplementary_reads + other.supplementary_reads,
            total_mapped=self.total_mapped + other.total_mapped,
            primary_mapped=self.primary_mapped + other.primary_mapped,
        )


@dataclass
class SampleQCResult:
    name: str
    group: str

    flagstat: FlagstatResult
    nanoplot: pd.DataFrame

    @property
    def label(self) -> str:
        return f"{self.group}/{self.name}"


@dataclass
class QCResult:
    sample_results: list[SampleQCResult]
    samples_df: pd.DataFrame
