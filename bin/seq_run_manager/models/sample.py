from dataclasses import dataclass
from pathlib import Path


@dataclass
class Sample:
    run_number: int
    replicate_number: int
    replicate_dir: Path
