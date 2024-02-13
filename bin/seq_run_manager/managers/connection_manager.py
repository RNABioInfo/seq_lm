from ..models.run_config import RunConfig
from typing import List

from pprint import pprint
import minknow_api as mk
import secrets
import string
from typing import Optional
import time


class ConnectionManager:
    manager: mk.manager.Manager

    def __init__(self, manager: mk.manager.Manager) -> None:
        self.manager = manager

    @classmethod
    def connected_with(cls, run_config: RunConfig) -> "ConnectionManager":
        manager = cls.__connect_to_minknow(run_config)
        return cls(manager)

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

    def __create_simulated_position(self, name: Optional[str] = None) -> str:
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

    def __get_simulated_position_by_name(
        self, name: str, retries: int = 3
    ) -> Optional[mk.manager.FlowCellPosition]:
        print(f"Looking for simulated position {name}")
        for _ in range(retries):
            for device in self.manager.flow_cell_positions():
                if device.name == name and device.is_simulated:
                    return device
            time.sleep(1)  # Wait for 1 second before retrying
        return None

    def remove_all_simulated_positions(self):
        positions = self.manager.flow_cell_positions()

        for position in positions:
            if position.is_simulated:
                self.manager.remove_simulated_device(str(position.name))

    def get_available_positions(self) -> List[mk.manager.FlowCellPosition]:
        return list(self.manager.flow_cell_positions())

    def print_available_positions(self):
        positions = self.manager.flow_cell_positions()

        for position in positions:
            pprint(position)

    def get_sequencing_positions(
        self, run_config: RunConfig
    ) -> List[mk.manager.FlowCellPosition]:
        if run_config.simulate_run:
            simulated_positions: List[mk.manager.FlowCellPosition] = []

            for _ in range(run_config.replicate_count):
                position_name = self.__create_simulated_position()
                simulated_position = self.__get_simulated_position_by_name(
                    position_name
                )

                if simulated_position is None:
                    raise Exception(
                        f"Could not retrieve simulated position {position_name}"
                    )

                simulated_positions.append(simulated_position)

            return simulated_positions

        if run_config.position_ids is not None and run_config.flow_cell_ids is not None:
            raise Exception(
                "You can only specify either a position_id or a flow_cell_id"
            )

        if run_config.position_ids is not None:
            if len(run_config.position_ids) != run_config.replicate_count:
                raise Exception(
                    f"Number of positions ({len(run_config.position_ids)}) does not match the number of replicates ({run_config.replicate_count})"
                )
            return [
                self.__get_sequencing_position_by_position_id_throws(id)
                for id in run_config.position_ids
            ]

        if run_config.flow_cell_ids is not None:
            if len(run_config.flow_cell_ids) != run_config.replicate_count:
                raise Exception(
                    f"Number of flow cells ({len(run_config.flow_cell_ids)}) does not match the number of replicates ({run_config.replicate_count})"
                )
            return [
                self.__get_sequencing_position_by_flow_cell_id_throws(id)
                for id in run_config.flow_cell_ids
            ]

        positions = list(self.manager.flow_cell_positions())

        if len(positions) != run_config.replicate_count:
            raise Exception(
                f"Number of available positions ({len(positions)}) does not match the number of replicates ({run_config.replicate_count})"
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

    def connect_to_positions(
        self, run_config: RunConfig, retries: int = 3
    ) -> List[mk.Connection]:
        positions = self.get_sequencing_positions(run_config)

        connections: List[mk.Connection] = []

        for position in positions:
            for _ in range(retries):
                try:
                    if not position.running:
                        raise Exception(f"Position {position.name} is not running.")

                    connections.append(position.connect())
                    break
                except Exception as e:
                    print(
                        f"Failed to connect to position {position.name}. Retrying... Info: {e}"
                    )
                    time.sleep(5)

        if len(connections) != run_config.replicate_count:
            raise Exception(
                f"Could not connect to all positions. Expected {run_config.replicate_count}, got {len(connections)}"
            )
        return connections
