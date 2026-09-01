from dataclasses import dataclass

from .connection_config import ConnectionConfig
from .manager_error import ManagerError


@dataclass(kw_only=True)
class StopAcquisitionConfig(ConnectionConfig):
    run_id: str

    def validate(self) -> None:
        super().validate()
        if not self.run_id.strip():
            raise ManagerError("Run ID is required")
