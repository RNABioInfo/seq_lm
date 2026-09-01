from dataclasses import dataclass
from pathlib import Path

from .connection_config import ConnectionConfig
from .manager_error import ManagerError
from .sample import Sample


@dataclass(kw_only=True)
class StartRunConfig(ConnectionConfig):
    experiment_id: str
    flow_cell_ids: list[str] | None
    position_ids: list[str] | None
    kit: str
    reference_genome_path: Path | None
    sampling_regions_path: Path | None
    adaptive_sampling_mode: str | None
    basecall_model: str | None
    min_qscore: float | None
    output_chunk_size: int
    samples: list[Sample]
    simulate_run: bool

    @property
    def sample_count(self) -> int:
        return len(self.samples)

    def validate(self) -> None:
        super().validate()

        if not self.experiment_id.strip():
            raise ManagerError("Experiment ID is required")
        if not self.kit.strip():
            raise ManagerError("Kit is required")

        if self.reference_genome_path:
            allowed_extensions = {".fasta", ".fa", ".fna", ".mmi"}
            if self.reference_genome_path.suffix.lower() not in allowed_extensions:
                raise ManagerError(
                    "Reference genome must be a .fasta, .fa, .fna or .mmi file"
                )

        if self.sampling_regions_path:
            if self.sampling_regions_path.suffix.lower() != ".bed":
                raise ManagerError("Sampling regions must be a .bed file")
            if not self.reference_genome_path:
                raise ManagerError("Sampling regions require a reference genome")

        if self.adaptive_sampling_mode and not self.sampling_regions_path:
            raise ManagerError("Adaptive sampling requires a sampling-regions BED file")
        if self.basecall_model and self.basecall_model.endswith(".cfg"):
            raise ManagerError(
                "Basecall model must be a Dorado simplex model name, not a legacy .cfg file"
            )
        if self.min_qscore is not None and self.min_qscore < 0:
            raise ManagerError("Minimum Q-score must be non-negative")
        if self.output_chunk_size < 1:
            raise ManagerError("Output chunk size must be greater than zero")
        if not self.samples:
            raise ManagerError("No samples provided")
