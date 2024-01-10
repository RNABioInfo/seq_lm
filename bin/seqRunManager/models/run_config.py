from typing import Optional


class RunConfig:
    host: str
    port: int
    certificate_path: str
    key_path: str
    flow_cell_id: Optional[str]
    position_id: Optional[str]
    sample_id: str
    experiment_id: str
    flow_cell_product_code: str
    kit: str
    reference_genome_path: str
    sampling_regions_path: Optional[str]
    adaptive_sampling_mode: Optional[str]
    basecall_config: str
    output_chunk_size: int
    output_dir: str

    def __init__(
        self,
        host: str,
        port: int,
        certificate_path: str,
        key_path: str,
        flow_cell_id: str,
        position_id: str,
        sample_id: str,
        experiment_id: str,
        flow_cell_product_code: str,
        kit: str,
        reference_genome_path: str,
        sampling_regions_path: Optional[str],
        adaptive_sampling_mode: Optional[str],
        basecall_config: str,
        output_chunk_size: int,
        output_dir: str,
    ):
        self.host = host
        self.port = port
        self.certificate_path = certificate_path
        self.key_path = key_path
        self.flow_cell_id = flow_cell_id
        self.position_id = position_id
        self.sample_id = sample_id
        self.experiment_id = experiment_id
        self.flow_cell_product_code = flow_cell_product_code
        self.kit = kit
        self.reference_genome_path = reference_genome_path
        self.sampling_regions_path = sampling_regions_path
        self.adaptive_sampling_mode = adaptive_sampling_mode
        self.basecall_config = basecall_config
        self.output_chunk_size = output_chunk_size
        self.output_dir = output_dir
