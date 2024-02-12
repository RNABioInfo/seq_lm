from ..models.acquisition import Acquisition
import time


class RunManager:
    active_acquisitions: list[Acquisition]

    def __init__(self, active_acquisitions: list[Acquisition]) -> None:
        self.active_acquisitions = active_acquisitions

    def watch_acquisitions_for_stop(self) -> None:
        while True:
            active_acquisitions = False

            for acquisition in self.active_acquisitions:
                if acquisition.is_stopped:
                    continue

                if acquisition.should_stop():
                    try:
                        acquisition.stop_run_throws()
                        continue
                    except Exception as e:
                        print(e)

                active_acquisitions = True

            if not active_acquisitions:
                break

            time.sleep(5)
