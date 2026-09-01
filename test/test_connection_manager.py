import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from seq_run_manager.managers.connection_manager import ConnectionManager


class ConnectionManagerCredentialsTest(unittest.TestCase):
    def test_client_certificates_disable_local_token_discovery(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            certificate = root / "client.pem"
            private_key = root / "client.key"
            ca_certificate = root / "ca.crt"
            certificate.write_bytes(b"client certificate")
            private_key.write_bytes(b"private key")
            ca_certificate.write_bytes(b"server CA")
            config = SimpleNamespace(
                host="127.0.0.1",
                port=9501,
                client_certificate_path=certificate,
                client_private_key_path=private_key,
                ca_certificate_path=ca_certificate,
            )

            with patch(
                "seq_run_manager.managers.connection_manager.mk_manager.Manager"
            ) as manager:
                ConnectionManager._ConnectionManager__connect_to_minknow(config)

            arguments = manager.call_args.kwargs
            self.assertEqual(
                arguments["client_certificate_chain"], certificate.read_bytes()
            )
            self.assertEqual(arguments["client_private_key"], private_key.read_bytes())
            self.assertEqual(arguments["ca_certificate"], ca_certificate.read_bytes())
            self.assertEqual(arguments["environ"]["MINKNOW_API_USE_LOCAL_TOKEN"], "0")


if __name__ == "__main__":
    unittest.main()
