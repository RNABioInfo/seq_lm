import time

import minknow_api as mk
from minknow_api.tools import protocols

from ..managers.connection_manager import ConnectionManager
from ..managers.sequencing_protocol_manager import SequencingProtocolManager
from ..models.acquisition import Acquisition
from ..models.manager_error import ManagerError
from ..models.run_config import RunConfig


class RunManager:
    connection_manager: ConnectionManager
    active_acquisitions: list[Acquisition]
    # For more information on the protocol states, see the minknow_api.proto.protocol.ProtocolState enum
    __permitted_protocol_states: list [int]

    def __init__(self, connection_manager: ConnectionManager) -> None:
        self.connection_manager = connection_manager
        self.active_acquisitions = []
        self.__permitted_protocol_states = [
            0,
            4,
            5,
            1,
        ]  # 0: PROTOCOL_RUNNING, 4: PROTOCOL_WAITING_FOR_TEMPERATURE, 5: PROTOCOL_WAITING_FOR_ACQUISITION, 1: PROTOCOL_COMPLETED

    def __get_product_code(self, connection: mk.Connection) -> str:
        flow_cell_info = connection.device.get_flow_cell_info() # type: ignore
        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code
        return product_code

    def __get_uniform_product_code(self, connections: list[mk.Connection]) -> str:
        product_code = self.__get_product_code(connections[0])
        for connection in connections:
            other_product_code = self.__get_product_code(connection)
            if other_product_code != product_code:
                raise ManagerError("All flow cells must have the same product code")
        return product_code

    def __start_acquisitions(self, run_config: RunConfig) -> list[Acquisition]:
        connections: list[mk.Connection] = self.connection_manager.connect_to_positions(
            run_config
        )

        product_code = self.__get_uniform_product_code(connections)

        protocol: mk.protocol_pb2.ProtocolInfo | None = SequencingProtocolManager.get_sequencing_protocol(  # type: ignore
            connections[0], product_code, run_config.kit
        )

        if protocol is None:
            available_protocols = SequencingProtocolManager.get_available_protocols(
                connections[0]
            )
            available_kits = []
            for protocol in available_protocols:
                tags = protocol.tags  # type: ignore
                if tags["flow cell"].string_value == product_code:
                    available_kits.append(tags["kit"].string_value)

            raise ManagerError(
                f"No protocol identifier found. \nAvailable kits for this flow cell: {available_kits}"
            )

        simplex_model, min_qscore = self.__resolve_basecalling_settings(
            run_config, connections[0], protocol, product_code
        )

        if len(connections) != len(run_config.samples):
            raise ManagerError(
                "Number of connected devices does not match number of requested samples"
            )

        acquisitions: list[Acquisition] = []

        for connection, sample in zip(connections, run_config.samples):
            seq_id = SequencingProtocolManager.start_sequencing_protocol(
                connection,
                protocol,
                sample.id,
                sample.replicate_dir,
                run_config,
                simplex_model,
                min_qscore,
            )
            acquisitions.append(
                Acquisition(sample=sample, id=seq_id, connection=connection)
            )
            print(f"Started sequencing with id {seq_id} for {sample.id}")

        return acquisitions

    def __resolve_basecalling_settings(
        self,
        run_config: RunConfig,
        position_connection: mk.Connection,
        protocol: mk.protocol_pb2.ProtocolInfo,  # type: ignore
        product_code: str,
    ) -> tuple[str, float]:
        sample_rate = protocol.tags["sample rate"].int_value  # type: ignore
        configurations = list(
            self.connection_manager.manager.find_basecall_configurations(
                flow_cell_product_code=product_code,
                sequencing_kit=run_config.kit,
                sampling_rate=sample_rate,
                include_outdated=False,
            )
        )
        if not configurations:
            raise ManagerError(
                "No current basecalling models are available for "
                f"flow cell {product_code}, kit {run_config.kit}, and "
                f"sampling rate {sample_rate}"
            )

        try:
            if run_config.basecall_model:
                simplex_model = protocols.find_simplex_model(
                    configurations, run_config.basecall_model
                )
            else:
                _, simplex_model = protocols.find_default_simplex_model(
                    position_connection,
                    run_config.kit,
                    sample_rate,
                    configurations,
                )
        except RuntimeError as error:
            raise ManagerError(str(error)) from error

        min_qscore = (
            run_config.min_qscore
            if run_config.min_qscore is not None
            else simplex_model.default_q_score_cutoff # type: ignore
        )
        return simplex_model.name, min_qscore # type: ignore

    def __watch_acquisitions_status(self) -> None:
        while True:
            active_acquisitions = False

            for acquisition in self.active_acquisitions:

                if acquisition.is_stopped:
                    continue

                # Check if the acquisition is still running. Retry until the protocol update contains a state.
                for event in SequencingProtocolManager.stream_current_protocol_updates(
                    acquisition.connection
                ):
                    if hasattr(event, "state"):
                        if not event.state in self.__permitted_protocol_states:
                            raise ManagerError(
                                f"Unexpected protocol state: {event.state}. Check the minknow log for more information."
                            )
                        break

                if acquisition.should_stop():
                    try:
                        acquisition.stop_run_throws()
                        continue
                    except Exception as e:  # noqa: BLE001
                        print(
                            f"Error stopping acquisition: {e}. Please stop run: {event.run_id} manually."
                        )

                active_acquisitions = True

            if not active_acquisitions:
                break

            time.sleep(5)

    def start_run_watcher(self, run_config: RunConfig) -> None:
        self.active_acquisitions = self.__start_acquisitions(run_config)
        self.__watch_acquisitions_status()

        if run_config.simulate_run:
            self.connection_manager.remove_all_simulated_positions()

    # def stop_run(self, )