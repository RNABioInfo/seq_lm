from dataclasses import dataclass
from pathlib import Path


@dataclass
class CertificateConfig:
    output_directory: Path
    ca_certificate_source: Path | None
    minknow_client_certs_directory: Path | None
    common_name: str
    valid_days: int
    key_size: int
    force: bool
