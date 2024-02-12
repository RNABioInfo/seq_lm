from ..models.run_config import RunConfig
from ..models.acquisition import Acquisition
from .sequencing_protocol_manager import SequencingProtocolManager

from pprint import pprint
import minknow_api as mk
import secrets
import string
from typing import Optional


class ConnectionManager:
    number: int = 5

    run_config: RunConfig
    manager: mk.manager.Manager

    def __init__(self, run_config: RunConfig):
        self.run_config = run_config
        self.manager = ConnectionManager.__connect_to_minknow(run_config)

    @staticmethod
    def __connect_to_minknow(run_config: RunConfig) -> mk.manager.Manager:
        certificate_bytes: Optional[bytes] = None
        if run_config.certificate_path:
            certificate_bytes = open(run_config.certificate_path, "rb").read()

        key_bytes: Optional[bytes] = None
        if run_config.key_path:
            key_bytes = open(run_config.key_path, "rb").read()

        return mk.manager.Manager(
            host=run_config.host or "127.0.0.1",
            port=run_config.port,
            client_certificate_chain=certificate_bytes,
            client_private_key=key_bytes,
        )

    def disconnect(self) -> None:
        self.manager.close()

    def create_simulated_position(self, name: Optional[str] = None) -> str:
        def generate_random_string(length=8) -> str:
            characters = string.ascii_letters + string.digits
            random_string = "".join(secrets.choice(characters) for _ in range(length))
            return random_string

        if name is None:
            name = generate_random_string()

        if self.__get_sequencing_position_by_position_id(name) is not None:
            return name

        self.manager.add_simulated_device(
            name,
            mk.manager_pb2.SimulatedDeviceType.SIMULATED_P2,  # type: ignore
        )

        return name

    def remove_all_simulated_positions(self):
        positions = self.manager.flow_cell_positions()

        for position in positions:
            if position.is_simulated:
                self.manager.remove_simulated_device(str(position.name))

    def get_available_positions(self) -> list[mk.manager.FlowCellPosition]:
        return list(self.manager.flow_cell_positions())

    def print_available_positions(self):
        positions = self.manager.flow_cell_positions()

        for position in positions:
            pprint(position)

    def get_sequencing_positions_for_config(
        self,
    ) -> list[mk.manager.FlowCellPosition]:
        if (
            self.run_config.position_ids is not None
            and self.run_config.flow_cell_ids is not None
        ):
            raise Exception(
                "You can only specify either a position_id or a flow_cell_id"
            )

        if self.run_config.position_ids is not None:
            if len(self.run_config.position_ids) != self.run_config.replicate_count:
                raise Exception(
                    f"Number of positions ({len(self.run_config.position_ids)}) does not match the number of replicates ({self.run_config.replicate_count})"
                )
            return [
                self.__get_sequencing_position_by_position_id_throws(id)
                for id in self.run_config.position_ids
            ]

        if self.run_config.flow_cell_ids is not None:
            if len(self.run_config.flow_cell_ids) != self.run_config.replicate_count:
                raise Exception(
                    f"Number of flow cells ({len(self.run_config.flow_cell_ids)}) does not match the number of replicates ({self.run_config.replicate_count})"
                )
            return [
                self.__get_sequencing_position_by_flow_cell_id_throws(id)
                for id in self.run_config.flow_cell_ids
            ]

        positions = list(self.manager.flow_cell_positions())

        if len(positions) != self.run_config.replicate_count:
            raise Exception(
                f"Number of available positions ({len(positions)}) does not match the number of replicates ({self.run_config.replicate_count})"
            )

        return positions

    def __get_sequencing_position_by_position_id_throws(
        self, position_id: str
    ) -> mk.manager.FlowCellPosition:
        positions = self.manager.flow_cell_positions()
        for position in positions:
            if position.name == position_id:
                return position

        raise Exception(f"Position with id {position_id} not found")

    def __get_sequencing_position_by_position_id(
        self, position_id: str
    ) -> Optional[mk.manager.FlowCellPosition]:
        positions = self.manager.flow_cell_positions()
        for position in positions:
            if position.name == position_id:
                return position

        return None

    def __get_sequencing_position_by_flow_cell_id_throws(
        self, flow_cell_id: str
    ) -> mk.manager.FlowCellPosition:
        positions = self.manager.flow_cell_positions()
        for position in positions:
            position_connection = position.connect()
            flow_cell_info = position_connection.device.get_flow_cell_info()  # type: ignore
            if (
                flow_cell_info.flow_cell_id == flow_cell_id
                or flow_cell_info.user_specified_flow_cell_id == flow_cell_id
            ):
                return position

        raise Exception(f"Position with id {flow_cell_id} not found")

    def __connect_to_positions_for_config(self) -> list[mk.Connection]:
        positions = self.get_sequencing_positions_for_config()

        connections: list[mk.Connection] = []

        for position in positions:
            if not position.running:
                raise Exception(f"Position {position.name} is not running")

            connections.append(position.connect())

        return connections

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

    def start_run(self) -> list[Acquisition]:
        connections: list[mk.Connection] = self.__connect_to_positions_for_config()
        product_code = self.__get_uniform_product_code(connections)

        protocol: Optional[mk.protocol_pb2.ProtocolInfo] = SequencingProtocolManager.get_sequencing_protocol(  # type: ignore
            connections[0], product_code, self.run_config.kit
        )

        if protocol is None:
            raise Exception("No protocol identifier found")

        if len(connections) != len(self.run_config.samples):
            raise Exception(
                "Number of connected devices does not match number of requested samples"
            )

        acquisitions: list[Acquisition] = []

        for connection, sample in zip(connections, self.run_config.samples):
            seq_id = SequencingProtocolManager.start_sequencing_protocol(
                connection, protocol, sample.id, sample.replicate_dir, self.run_config
            )
            acquisitions.append(
                Acquisition(sample=sample, id=seq_id, connection=connection)
            )
            print(f"Started sequencing with id {seq_id} for {sample.id}")

        return acquisitions
