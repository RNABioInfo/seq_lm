from .sample import Sample
from datetime import datetime
import minknow_api as mk
from typing import Optional


class Acquisition:
    sample: Sample
    id: str
    connection: mk.Connection
    start_date: datetime
    end_date: Optional[datetime]
    is_stopped: bool = False

    def __init__(self, sample: Sample, id: str, connection: mk.Connection) -> None:
        self.sample = sample
        self.id = id
        self.connection = connection
        self.start_date = datetime.now()

    def should_stop(self) -> bool:
        stop_file_path = self.sample.replicate_dir / "STOP"

        if stop_file_path.exists():
            return True

        return False

    def stop_run_throws(self) -> None:
        try:
            self.connection.protocol.stop_protocol()  # type: ignore
            self.end_date = datetime.now()
            print(f"Run {self.sample.id} stopped")
        except Exception as e:
            raise Exception(f"Could not stop the run. Reason: {e}")
