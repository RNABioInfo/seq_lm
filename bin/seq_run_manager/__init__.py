#!/bin/python3

import socket
from collections.abc import Sequence

from .managers.certificate_manager import CertificateManager, CertificateSetupError
from .models.certificate_config import CertificateConfig
from .models.manager_error import ManagerError
from .models.start_run_config import StartRunConfig
from .models.stop_acquisition_config import StopAcquisitionConfig
from .utils.argument_parser import ArgumentParser

__version__ = "0.0.1"
_package_name = "seq_run_manager"


def main(arguments: Sequence[str] | None = None) -> None:
    config = ArgumentParser.parse_cli_arguments(arguments)
    print(f"Running {_package_name} version {__version__}")

    if isinstance(config, CertificateConfig):
        setup_certificates(config)
        return

    from .managers.connection_manager import ConnectionManager
    from .managers.run_manager import RunManager

    try:
        check_server(config.host, config.port)
    except ManagerError as error:
        raise SystemExit(f"seq-run-manager failed: {error}") from error
    
    connection_manager = ConnectionManager.connected_with(config)

    try:
        run_manager = RunManager(connection_manager)
        if isinstance(config, StartRunConfig):
            if config.simulate_run:
                connection_manager.remove_all_simulated_positions()
            run_manager.start_run_watcher(config)

        elif isinstance(config, StopAcquisitionConfig):
            run_manager.stop_run(config.run_id)
            print(f"Stopped acquisition with run ID {config.run_id}")
            
    except ManagerError as error:
        raise SystemExit(f"seq-run-manager failed: {error}") from error
    finally:
        try:
            if isinstance(config, StartRunConfig) and config.simulate_run:
                connection_manager.remove_all_simulated_positions()
        finally:
            connection_manager.disconnect()


def setup_certificates(config: CertificateConfig) -> None:
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

    print("Connection arguments for 'seq-run-manager start' or 'stop':")
    print(f"  --client-certificate-path {client_certificate}")
    print(f"  --client-private-key-path {client_private_key}")
    print(f"  --ca-certificate-path {ca_certificate}")


def check_server(address: str, port: int) -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as connection:
        try:
            connection.connect((address, port))
            print(f"Server is running at {address}:{port}")
        except OSError as error:
            raise ManagerError(
                f"Server not reachable at {address}:{port}. Error: {error}"
            ) from error
