from dataclasses import dataclass
from pathlib import Path


@dataclass
class Sample:
    run_number: int
    replicate_number: int
    replicate_dir: Path

    @property
    def id(self):
        return f"run_{self.run_number}_replicate_{self.replicate_number}"
