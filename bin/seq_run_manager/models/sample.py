from dataclasses import dataclass
from pathlib import Path


@dataclass
class Sample:
    alias: str
    group: str
    bam_dir: Path
    is_live: bool = True
    order: int | None = None

    @property
    def id(self) -> str:
        return self.alias
