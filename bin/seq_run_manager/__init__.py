import socket
from .models.run_config import RunConfig
from .utils.argument_parser import ArgumentParser
from .managers.connection_manager import ConnectionManager
from .managers.run_manager import RunManager
import sys

__version__ = "0.0.1"
_package_name = "seq_run_manager"


def main():
    print(f"Running {_package_name} version {__version__}")

    run_config: RunConfig = ArgumentParser.parse_cli_arguments()

    check_server(run_config.host, run_config.port)

    connection_manager: ConnectionManager = ConnectionManager.connected_with(run_config)
    connection_manager.print_available_positions()
    connection_manager.remove_all_simulated_positions()

    run_manager: RunManager = RunManager(connection_manager)

    if run_config.simulate_run:
        sys.exit(0)

    run_manager.start_run(run_config)

    connection_manager.disconnect()


def check_server(address, port):
    # Create a TCP socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.connect((address, port))
            print(f"Server is running at {address}:{port}")
        except socket.error as e:
            raise Exception(f"Server not reachable at {address}:{port}. Error: {e}")
