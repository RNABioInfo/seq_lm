from datetime import datetime

import minknow_api as mk

from .manager_error import ManagerError
from .sample import Sample


class Acquisition:
    sample: Sample
    id: str
    connection: mk.Connection
    start_date: datetime
    end_date: datetime | None
    is_stopped: bool = False

    def __init__(self, sample: Sample, id: str, connection: mk.Connection) -> None:
        self.sample = sample
        self.id = id
        self.connection = connection
        self.start_date = datetime.now()  # noqa: DTZ005

    def should_stop(self) -> bool:
        stop_file_path = self.sample.replicate_dir / "STOP"

        return bool(stop_file_path.exists())

    def stop_run_throws(self) -> None:
        try:
            self.connection.protocol.stop_protocol()  # type: ignore
            self.end_date = datetime.now()  # noqa: DTZ005
            print(f"Run {self.sample.id} stopped")
        except Exception as e:  # noqa: BLE001
            raise ManagerError(f"Could not stop the run. Reason: {e}")
