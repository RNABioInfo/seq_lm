import socket
from models.run_config import RunConfig
from utils.argument_parser import ArgumentParser
from managers.connection_manager import ConnectionManager
from managers.sequencing_protocol_manager import SequencingProtocolManager

# "/Applications/MinKNOW.app/Contents/Resources/conf/rpc-client-certs/minknow_cert.pem"
# "/Users/christopherphd/Documents/projects/bios/seqLM/minknow_key.pem"


def main(args=None):
    run_config: RunConfig = ArgumentParser.parse_cli_arguments()
    connection_manager: ConnectionManager = ConnectionManager(run_config)

    # connection_manager.get_sequencing_position_for_config()
    connection_manager.create_simulated_position("simulated")

    # wait for two seconds
    import time

    time.sleep(2)
    position = connection_manager.get_sequencing_position_for_config()
    connection = position.connect()

    # protocol = SequencingProtocolManager.get_sequencing_protocol(
    #     connection, "FLO-PRO004RA", "SQK-RNA004", "rna_rp4_130bps_hac_prom.cfg"
    # )

    run_id: str = SequencingProtocolManager.start_sequencing_protocol(
        connection, run_config
    )
    run_info = SequencingProtocolManager.stream_acquisition_output(connection, run_id)

    time.sleep(30)

    connection_manager.remove_all_simulated_positions()


def check_server(address, port):
    # Create a TCP socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.connect((address, port))
            print(f"Server is running at {address}:{port}")
        except socket.error as e:
            raise Exception(f"Server not reachable at {address}:{port}. Error: {e}")


if "__main__" == __name__:
    check_server("127.0.0.1", 9501)
    main()
