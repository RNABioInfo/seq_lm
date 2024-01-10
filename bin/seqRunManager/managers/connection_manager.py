from models.run_config import RunConfig

from pprint import pprint
import minknow_api as mk
import secrets
import string
from typing import Optional


class ConnectionManager:
    run_config: RunConfig
    manager: mk.manager.Manager

    def __init__(self, run_config: RunConfig):
        self.run_config = run_config
        self.manager = ConnectionManager.__connect_to_minknow(run_config)

    @staticmethod
    def __connect_to_minknow(run_config: RunConfig) -> mk.manager.Manager:
        certificate_bytes: bytes = open(run_config.certificate_path, "rb").read()
        key_bytes: bytes = open(run_config.key_path, "rb").read()

        return mk.manager.Manager(
            host=run_config.host,
            port=run_config.port,
            client_certificate_chain=certificate_bytes,
            client_private_key=key_bytes,
        )

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

    def print_available_positions(self):
        positions = self.manager.flow_cell_positions()

        for position in positions:
            pprint(position)

    def get_sequencing_position_for_config(self) -> mk.manager.FlowCellPosition:
        if (
            self.run_config.position_id is not None
            and self.run_config.flow_cell_id is not None
        ):
            raise Exception(
                "You can only specify either a position_id or a flow_cell_id"
            )

        if self.run_config.position_id is not None:
            return self.__get_sequencing_position_by_position_id_throws(
                self.run_config.position_id
            )

        if self.run_config.flow_cell_id is not None:
            return self.__get_sequencing_position_by_flow_cell_id(
                self.run_config.flow_cell_id
            )

        raise Exception(
            f"Could not find a sequencing position for the given config: {self.run_config}"
        )

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

    def __get_sequencing_position_by_flow_cell_id(
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

    def connect_to_position_by_config(self) -> mk.Connection:
        position = self.get_sequencing_position_for_config()

        if not position.running:
            raise Exception(f"Position {position.name} is not running")

        return position.connect()
