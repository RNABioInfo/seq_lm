from typing import Optional, List
from dataclasses import dataclass
from .sample import Sample


@dataclass
class RunConfig:
    host: Optional[str]
    port: Optional[int]
    certificate_path: Optional[str]
    key_path: Optional[str]
    flow_cell_ids: Optional[List[str]]
    position_ids: Optional[List[str]]
    experiment_id: str
    run_number: int
    replicate_count: int
    kit: str
    reference_genome_path: str
    sampling_regions_path: Optional[str]
    adaptive_sampling_mode: Optional[str]
    basecall_config: str
    output_chunk_size: int
    samples: List[Sample]
    simulate_run: bool
