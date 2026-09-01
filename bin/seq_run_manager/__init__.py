#!/bin/python3

import socket
import sys

from .managers.certificate_manager import CertificateManager, CertificateSetupError
from .managers.sequencing_protocol_manager import SequencingProtocolManager
from .models.manager_error import ManagerError
from .utils.argument_parser import ArgumentParser

__version__ = "0.0.1"
_package_name = "seq_run_manager"


def main():
    print(f"Running {_package_name} version {__version__}")

    if len(sys.argv) > 1 and sys.argv[1] == "setup-certificates":
        setup_certificates()
        return

    from .managers.connection_manager import ConnectionManager
    from .managers.run_manager import RunManager
    from .models.run_config import RunConfig

    run_config: RunConfig = ArgumentParser.parse_cli_arguments()

    check_server(run_config.host, run_config.port)

    connection_manager: ConnectionManager = ConnectionManager.connected_with(run_config)
    connection_manager.print_available_positions()
    connection_manager.remove_all_simulated_positions()

    positions = connection_manager.get_available_positions()

    for pos in positions:
        connection = pos.connect()
        active = SequencingProtocolManager.get_currently_active_protocol(connection)
        print(f"active: {active}")

    return

    run_manager: RunManager = RunManager(connection_manager)

    if run_config.simulate_run:
        sys.exit(0)

    run_manager.start_run_watcher(run_config)

    connection_manager.disconnect()


def setup_certificates():
    config = ArgumentParser.parse_certificate_setup_arguments(sys.argv[2:])
    try:
        (
            client_certificate,
            client_private_key,
            client_ca_certificate,
            ca_certificate,
        ) = CertificateManager.setup(config)
    except CertificateSetupError as error:
        raise SystemExit(f"Certificate setup failed: {error}") from error

    print("MinKNOW credentials created and validated:")
    print(f"  Client certificate: {client_certificate}")
    print(f"  Client private key: {client_private_key}")
    print(f"  Client CA certificate: {client_ca_certificate}")
    print(f"  CA certificate: {ca_certificate}")
    print(
        "  Client certificate SHA-256: "
        f"{CertificateManager.certificate_fingerprint(client_certificate)}"
    )
    if config.minknow_client_certs_directory is not None:
        installed_certificate = (
            config.minknow_client_certs_directory.expanduser().resolve()
            / CertificateManager.INSTALLED_CLIENT_CERTIFICATE_NAME
        )
        print(f"  Installed public certificate: {installed_certificate}")
        print("Restart MinKNOW if it does not reload trusted client certificates.")
    else:
        print(
            "The public client certificate was not installed in MinKNOW. Run this "
            "command again with --minknow-client-certs-directory pointing to "
            "MinKNOW's conf/rpc-client-certs directory."
        )

    print("Connection arguments:")
    print(f"  --client_certificate_path {client_certificate}")
    print(f"  --client_private_key_path {client_private_key}")
    print(f"  --ca_certificate_path {ca_certificate}")


def check_server(address, port):
    # Create a TCP socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.connect((address, port))
            print(f"Server is running at {address}:{port}")
        except OSError as e:
            raise ManagerError(f"Server not reachable at {address}:{port}. Error: {e}")
