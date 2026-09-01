from dataclasses import dataclass
from pathlib import Path

from .manager_error import ManagerError


@dataclass(kw_only=True)
class ConnectionConfig:
    host: str
    port: int
    client_certificate_path: Path
    client_private_key_path: Path
    ca_certificate_path: Path

    def validate(self) -> None:
        if not self.host.strip():
            raise ManagerError("Host must not be empty")
        if not 1 <= self.port <= 65535:
            raise ManagerError("Port must be between 1 and 65535")

        credential_paths = (
            (self.client_certificate_path, "Client certificate"),
            (self.client_private_key_path, "Client private key"),
            (self.ca_certificate_path, "CA certificate"),
        )
        for path, description in credential_paths:
            if not path.is_file():
                raise ManagerError(f"{description} file does not exist: {path}")
