import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from seq_run_manager.managers.certificate_manager import (
    CertificateManager,
    CertificateSetupError,
)
from seq_run_manager.models.certificate_config import CertificateConfig


@unittest.skipUnless(shutil.which("openssl"), "OpenSSL is required")
class CertificateManagerTest(unittest.TestCase):
    def test_creates_valid_credentials_and_installs_only_public_certificate(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source_ca = root / "ca.crt"
            source_ca_key = root / "ca.key"
            output_directory = root / "credentials"
            trust_directory = root / "rpc-client-certs"
            trust_directory.mkdir()
            self._create_test_ca(source_ca, source_ca_key)

            config = CertificateConfig(
                output_directory=output_directory,
                ca_certificate_source=source_ca,
                minknow_client_certs_directory=trust_directory,
                common_name="certificate-manager-test",
                valid_days=1,
                key_size=2048,
                force=False,
            )
            certificate, private_key, client_ca, copied_ca = CertificateManager.setup(
                config
            )

            self.assertTrue(certificate.is_file())
            self.assertTrue(private_key.is_file())
            self.assertEqual(source_ca.read_bytes(), copied_ca.read_bytes())
            self.assertEqual(
                client_ca.read_bytes(),
                (trust_directory / "seq-run-manager.pem").read_bytes(),
            )
            self.assertEqual(stat.S_IMODE(private_key.stat().st_mode), 0o600)
            self.assertFalse((trust_directory / private_key.name).exists())

            with self.assertRaises(CertificateSetupError):
                CertificateManager.setup(config)

    def test_uses_windows_elevation_when_wsl_mount_rejects_direct_copy(self):
        certificate = Path("/tmp/minknow_cert.pem")
        installation_directory = Path("/mnt/c/Program Files/MinKNOW/conf")

        with (
            patch.object(shutil, "copyfile", side_effect=PermissionError),
            patch.object(CertificateManager, "_is_wsl", return_value=True),
            patch.object(
                CertificateManager,
                "_is_windows_mounted_path",
                return_value=True,
            ),
            patch.object(
                CertificateManager,
                "_install_client_certificate_with_windows_elevation",
            ) as elevated_install,
            patch.object(Path, "unlink"),
        ):
            CertificateManager._install_client_certificate(
                certificate, installation_directory, force=False
            )

        elevated_install.assert_called_once_with(
            certificate,
            installation_directory.resolve() / "seq-run-manager.pem",
            False,
        )

    @staticmethod
    def _create_test_ca(certificate: Path, private_key: Path) -> None:
        subprocess.run(
            (
                shutil.which("openssl"),
                "req",
                "-x509",
                "-newkey",
                "rsa:2048",
                "-nodes",
                "-days",
                "1",
                "-keyout",
                str(private_key),
                "-out",
                str(certificate),
                "-subj",
                "/CN=certificate-manager-test-ca",
            ),
            check=True,
            capture_output=True,
        )


if __name__ == "__main__":
    unittest.main()
