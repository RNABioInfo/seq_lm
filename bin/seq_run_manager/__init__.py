import socket
import time
from .models.run_config import RunConfig
from .utils.argument_parser import ArgumentParser
from .managers.connection_manager import ConnectionManager

__version__ = "0.0.1"
_package_name = "seq_run_manager"


def main():
    print(f"Running {_package_name} version {__version__}")

    run_config: RunConfig = ArgumentParser.parse_cli_arguments()

    check_server(run_config.host, run_config.port)

    connection_manager: ConnectionManager = ConnectionManager(run_config)
    connection_manager.remove_all_simulated_positions()

    connection_manager.print_available_positions()

    exit(0)

    connection_manager.start_run()

    time.sleep(20)

    connection_manager.remove_all_simulated_positions()


def check_server(address, port):
    # Create a TCP socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.connect((address, port))
            print(f"Server is running at {address}:{port}")
        except socket.error as e:
            raise Exception(f"Server not reachable at {address}:{port}. Error: {e}")
