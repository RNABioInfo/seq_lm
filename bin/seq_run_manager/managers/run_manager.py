from ..models.acquisition import Acquisition
from ..models.run_config import RunConfig
from ..managers.sequencing_protocol_manager import SequencingProtocolManager
from ..managers.connection_manager import ConnectionManager
import time
import minknow_api as mk


class RunManager:
    connection_manager: ConnectionManager
    active_acquisitions: list[Acquisition] = []

    def __init__(self, connection_manager: ConnectionManager) -> None:
        self.connection_manager = connection_manager

    def __get_product_code(self, connection: mk.Connection) -> str:
        flow_cell_info = connection.device.get_flow_cell_info()  # type: ignore
        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code
        return product_code

    def __get_uniform_product_code(self, connections: list[mk.Connection]) -> str:
        product_code = self.__get_product_code(connections[0])
        for connection in connections:
            other_product_code = self.__get_product_code(connection)
            if other_product_code != product_code:
                raise Exception("All flow cells must have the same product code")
        return product_code

    def __start_acquisitions(self, run_config: RunConfig) -> list[Acquisition]:
        connections: list[mk.Connection] = self.connection_manager.connect_to_positions(
            run_config
        )
        product_code = self.__get_uniform_product_code(connections)

        protocol: Optional[mk.protocol_pb2.ProtocolInfo] = SequencingProtocolManager.get_sequencing_protocol(  # type: ignore
            connections[0], product_code, run_config.kit
        )

        if protocol is None:
            raise Exception("No protocol identifier found")

        if len(connections) != len(run_config.samples):
            raise Exception(
                "Number of connected devices does not match number of requested samples"
            )

        acquisitions: list[Acquisition] = []

        for connection, sample in zip(connections, run_config.samples):
            seq_id = SequencingProtocolManager.start_sequencing_protocol(
                connection, protocol, sample.id, sample.replicate_dir, run_config
            )
            acquisitions.append(
                Acquisition(sample=sample, id=seq_id, connection=connection)
            )
            print(f"Started sequencing with id {seq_id} for {sample.id}")

        return acquisitions

    def __watch_acquisitions_for_stop(self) -> None:
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

    def start_run(self, run_config: RunConfig) -> None:
        self.active_acquisitions = self.__start_acquisitions(run_config)
        self.__watch_acquisitions_for_stop()
