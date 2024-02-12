from typing import Optional
from dataclasses import dataclass


@dataclass
class RunConfig:
    host: Optional[str]
    port: Optional[int]
    certificate_path: Optional[str]
    key_path: Optional[str]
    flow_cell_ids: Optional[list[str]]
    position_ids: Optional[list[str]]
    experiment_id: str
    run_number: int
    replicate_count: int
    kit: str
    reference_genome_path: str
    sampling_regions_path: Optional[str]
    adaptive_sampling_mode: Optional[str]
    basecall_config: str
    output_chunk_size: int
    output_dir: str
