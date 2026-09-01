from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .sample import Sample


@dataclass
class RunConfig:
    host: Optional[str]
    port: Optional[int]
    client_certificate_path: str
    client_private_key_path: str
    ca_certificate_path: str
    flow_cell_ids: Optional[List[str]]
    position_ids: Optional[List[str]]
    experiment_id: str
    run_id: str
    kit: str
    reference_genome_path: str
    sampling_regions_path: Optional[str]
    adaptive_sampling_mode: Optional[str]
    basecall_config: str
    output_chunk_size: int
    samples: List[Sample]
    simulate_run: bool

    @property
    def sample_count(self) -> int:
        return len(self.samples)

    def validate(self):
        if not Path(self.client_certificate_path).exists():
            raise Exception("Client certificate file does not exist")

        if not Path(self.client_private_key_path).exists():
            raise Exception("Client private key file does not exist")

        if not Path(self.ca_certificate_path).exists():
            raise Exception("CA certificate file does not exist")

        if self.experiment_id is None:
            raise Exception("Experiment ID is required")

        if self.run_id is None:
            raise Exception("Run ID is required")

        if self.kit is None:
            raise Exception("Kit is required")

        if self.reference_genome_path:
            allowed_reference_genome_extensions = [".fasta", ".fa", ".fna", ".mmi"]
            reference_genome_extension = Path(self.reference_genome_path).suffix
            if reference_genome_extension not in allowed_reference_genome_extensions:
                raise Exception(
                    "Reference genome must be a .fasta, .fa, .fna or .mmi file"
                )

        if Path(self.basecall_config).suffix != ".cfg":
            raise Exception("Basecall config must be a .cfg file")

        if len(self.samples) < 1:
            raise Exception("No samples provided")
