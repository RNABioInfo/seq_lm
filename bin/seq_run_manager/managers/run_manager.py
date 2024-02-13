from ..models.acquisition import Acquisition
from ..models.run_config import RunConfig
from ..managers.sequencing_protocol_manager import SequencingProtocolManager
from ..managers.connection_manager import ConnectionManager
from typing import List
import time
import minknow_api as mk
from pathlib import Path


class RunManager:
    connection_manager: ConnectionManager
    active_acquisitions: List[Acquisition] = []

    def __init__(self, connection_manager: ConnectionManager) -> None:
        self.connection_manager = connection_manager

    def __get_product_code(self, connection: mk.Connection) -> str:
        flow_cell_info = connection.device.get_flow_cell_info()  # type: ignore
        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code
        return product_code

    def __get_uniform_product_code(self, connections: List[mk.Connection]) -> str:
        product_code = self.__get_product_code(connections[0])
        for connection in connections:
            other_product_code = self.__get_product_code(connection)
            if other_product_code != product_code:
                raise Exception("All flow cells must have the same product code")
        return product_code

    def __start_acquisitions(self, run_config: RunConfig) -> List[Acquisition]:
        connections: List[mk.Connection] = self.connection_manager.connect_to_positions(
            run_config
        )

        product_code = self.__get_uniform_product_code(connections)

        self.__check_bascalling_config(run_config, product_code)

        protocol: Optional[mk.protocol_pb2.ProtocolInfo] = SequencingProtocolManager.get_sequencing_protocol(  # type: ignore
            connections[0], product_code, run_config.kit
        )

        if protocol is None:
            available_protocols = SequencingProtocolManager.get_available_protocols(
                connections[0]
            )
            available_kits = []
            for protocol in available_protocols:
                tags = protocol.tags
                if tags["flow cell"].string_value == product_code:
                    available_kits.append(tags["kit"].string_value)

            raise Exception(
                f"No protocol identifier found. \nAvailable kits for this flow cell: {available_kits}"
            )

        if len(connections) != len(run_config.samples):
            raise Exception(
                "Number of connected devices does not match number of requested samples"
            )

        acquisitions: List[Acquisition] = []

        for connection, sample in zip(connections, run_config.samples):
            seq_id = SequencingProtocolManager.start_sequencing_protocol(
                connection, protocol, sample.id, sample.replicate_dir, run_config
            )
            acquisitions.append(
                Acquisition(sample=sample, id=seq_id, connection=connection)
            )
            print(f"Started sequencing with id {seq_id} for {sample.id}")

        return acquisitions

    def __check_bascalling_config(
        self, run_config: RunConfig, uniform_product_code: str
    ) -> None:
        configs_by_flow_cell = mk.manager.Manager.basecaller(self.connection_manager.manager).rpc.list_configs_by_kit()  # type: ignore
        configs_for_run: list = (
            configs_by_flow_cell.flow_cell_configs[uniform_product_code]
            .kit_configs[run_config.kit]
            .configs
        )
        requested_config = str(Path(run_config.basecall_config).stem)
        if requested_config not in configs_for_run:
            raise Exception(
                f"Basecalling config {run_config.basecall_config} not available for flow cell {uniform_product_code} and kit {run_config.kit}. Available configs: {configs_for_run}"
            )

    def __watch_acquisitions_for_stop(self) -> None:
        while True:
            active_acquisitions = False

            for acquisition in self.active_acquisitions:
                for event in SequencingProtocolManager.stream_current_protocol_updates(
                    acquisition.connection
                ):
                    if hasattr(event, "state"):
                        print(event.state)
                    break

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

        if run_config.simulate_run:
            self.connection_manager.remove_all_simulated_positions()
