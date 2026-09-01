import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from ..models.certificate_setup_config import CertificateSetupConfig


class CertificateSetupError(Exception):
    pass


class CertificateManager:
    CLIENT_CERTIFICATE_NAME = "minknow_cert.pem"
    CLIENT_PRIVATE_KEY_NAME = "minknow_key.pem"
    CA_CERTIFICATE_NAME = "minknow_cert.crt"
    INSTALLED_CLIENT_CERTIFICATE_NAME = "seq-run-manager.pem"

    @classmethod
    def setup(cls, config: CertificateSetupConfig) -> tuple[Path, Path, Path]:
        openssl = shutil.which("openssl")
        if openssl is None:
            raise CertificateSetupError(
                "OpenSSL is required to create MinKNOW client credentials"
            )

        output_directory = config.output_directory.expanduser().resolve()
        output_directory.mkdir(mode=0o700, parents=True, exist_ok=True)

        client_certificate = output_directory / cls.CLIENT_CERTIFICATE_NAME
        client_private_key = output_directory / cls.CLIENT_PRIVATE_KEY_NAME
        ca_certificate = output_directory / cls.CA_CERTIFICATE_NAME
        targets = (client_certificate, client_private_key, ca_certificate)
        cls._ensure_replace_is_allowed(targets, config.force)
        if config.minknow_client_certs_directory is not None:
            cls._ensure_install_is_allowed(
                config.minknow_client_certs_directory, config.force
            )

        ca_source = cls._find_ca_certificate(config.ca_certificate_source)
        cls._validate_certificate(openssl, ca_source, "MinKNOW CA certificate")

        with tempfile.TemporaryDirectory(
            prefix=".certificate-setup-", dir=output_directory
        ) as temporary_directory:
            temporary_path = Path(temporary_directory)
            generated_certificate = temporary_path / cls.CLIENT_CERTIFICATE_NAME
            generated_private_key = temporary_path / cls.CLIENT_PRIVATE_KEY_NAME
            copied_ca_certificate = temporary_path / cls.CA_CERTIFICATE_NAME

            cls._generate_client_credentials(
                openssl=openssl,
                certificate=generated_certificate,
                private_key=generated_private_key,
                temporary_directory=temporary_path,
                common_name=config.common_name,
                valid_days=config.valid_days,
                key_size=config.key_size,
            )
            shutil.copyfile(ca_source, copied_ca_certificate)
            cls._validate_credentials(
                openssl,
                generated_certificate,
                generated_private_key,
                copied_ca_certificate,
            )

            os.replace(generated_certificate, client_certificate)
            os.replace(generated_private_key, client_private_key)
            os.replace(copied_ca_certificate, ca_certificate)

        client_certificate.chmod(0o644)
        client_private_key.chmod(0o600)
        ca_certificate.chmod(0o644)

        if config.minknow_client_certs_directory is not None:
            cls._install_client_certificate(
                client_certificate,
                config.minknow_client_certs_directory,
                config.force,
            )

        return client_certificate, client_private_key, ca_certificate

    @classmethod
    def certificate_fingerprint(cls, certificate: Path) -> str:
        openssl = shutil.which("openssl")
        if openssl is None:
            raise CertificateSetupError("OpenSSL is required to inspect certificates")
        result = cls._run(
            openssl,
            "x509",
            "-in",
            str(certificate),
            "-noout",
            "-fingerprint",
            "-sha256",
        )
        return result.stdout.strip().removeprefix("sha256 Fingerprint=")

    @classmethod
    def _find_ca_certificate(cls, configured_source: Path | None) -> Path:
        candidates: list[Path] = []
        if configured_source is not None:
            candidates.append(configured_source.expanduser())

        environment_source = os.environ.get("MINKNOW_TRUSTED_CA")
        if environment_source:
            candidates.append(Path(environment_source).expanduser())

        candidates.extend(
            (
                Path("/mnt/c/data/rpc-certs/minknow/ca.crt"),
                Path("/data/rpc-certs/minknow/ca.crt"),
                Path("/var/lib/minknow/data/rpc-certs/minknow/ca.crt"),
                Path("/Library/MinKNOW/data/rpc-certs/minknow/ca.crt"),
            )
        )

        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()

        searched = ", ".join(str(path) for path in candidates)
        raise CertificateSetupError(
            "Could not find the MinKNOW CA certificate. Pass "
            f"--ca-certificate-source explicitly. Searched: {searched}"
        )

    @staticmethod
    def _ensure_replace_is_allowed(targets: tuple[Path, ...], force: bool) -> None:
        existing = [path for path in targets if path.exists()]
        if existing and not force:
            paths = ", ".join(str(path) for path in existing)
            raise CertificateSetupError(
                f"Refusing to replace existing credentials: {paths}. Use --force."
            )

    @classmethod
    def _generate_client_credentials(
        cls,
        openssl: str,
        certificate: Path,
        private_key: Path,
        temporary_directory: Path,
        common_name: str,
        valid_days: int,
        key_size: int,
    ) -> None:
        if not common_name or any(
            character in common_name for character in ("/", "\n")
        ):
            raise CertificateSetupError(
                "The certificate Common Name must be non-empty and cannot contain '/' "
                "or a newline"
            )

        openssl_configuration = temporary_directory / "openssl-client.cnf"
        openssl_configuration.write_text(
            """[req]
distinguished_name = subject
x509_extensions = client_extensions
prompt = no

[subject]
CN = placeholder

[client_extensions]
basicConstraints = critical,CA:FALSE
keyUsage = critical,digitalSignature,keyEncipherment
extendedKeyUsage = clientAuth
"""
        )
        cls._run(
            openssl,
            "req",
            "-x509",
            "-newkey",
            f"rsa:{key_size}",
            "-sha256",
            "-nodes",
            "-days",
            str(valid_days),
            "-keyout",
            str(private_key),
            "-out",
            str(certificate),
            "-config",
            str(openssl_configuration),
            "-subj",
            f"/CN={common_name}",
        )

    @classmethod
    def _validate_credentials(
        cls,
        openssl: str,
        certificate: Path,
        private_key: Path,
        ca_certificate: Path,
    ) -> None:
        cls._validate_certificate(openssl, certificate, "client certificate")
        cls._validate_certificate(openssl, ca_certificate, "MinKNOW CA certificate")
        cls._run(openssl, "pkey", "-in", str(private_key), "-check", "-noout")

        purposes = cls._run(
            openssl, "x509", "-in", str(certificate), "-purpose", "-noout"
        ).stdout
        if "SSL client : Yes" not in purposes:
            raise CertificateSetupError(
                "Generated certificate is not valid for TLS client authentication"
            )

        certificate_key = cls._run(
            openssl, "x509", "-in", str(certificate), "-pubkey", "-noout"
        ).stdout
        private_key_public = cls._run(
            openssl, "pkey", "-in", str(private_key), "-pubout"
        ).stdout
        if certificate_key != private_key_public:
            raise CertificateSetupError(
                "Generated client certificate and private key do not match"
            )

    @classmethod
    def _validate_certificate(
        cls, openssl: str, certificate: Path, description: str
    ) -> None:
        try:
            cls._run(openssl, "x509", "-in", str(certificate), "-noout")
            cls._run(
                openssl,
                "x509",
                "-in",
                str(certificate),
                "-checkend",
                "0",
                "-noout",
            )
        except CertificateSetupError as error:
            raise CertificateSetupError(
                f"Invalid {description} at {certificate}: {error}"
            ) from error

    @classmethod
    def _install_client_certificate(
        cls, certificate: Path, configured_directory: Path, force: bool
    ) -> None:
        installation_directory = configured_directory.expanduser().resolve()
        installed_certificate = (
            installation_directory / cls.INSTALLED_CLIENT_CERTIFICATE_NAME
        )

        temporary_certificate = installation_directory / (
            f".{cls.INSTALLED_CLIENT_CERTIFICATE_NAME}.tmp"
        )
        try:
            shutil.copyfile(certificate, temporary_certificate)
            os.replace(temporary_certificate, installed_certificate)
            installed_certificate.chmod(0o644)
        finally:
            temporary_certificate.unlink(missing_ok=True)

    @classmethod
    def _ensure_install_is_allowed(
        cls, configured_directory: Path, force: bool
    ) -> None:
        installation_directory = configured_directory.expanduser().resolve()
        if not installation_directory.is_dir():
            raise CertificateSetupError(
                "MinKNOW client certificate directory does not exist: "
                f"{installation_directory}"
            )

        installed_certificate = (
            installation_directory / cls.INSTALLED_CLIENT_CERTIFICATE_NAME
        )
        if installed_certificate.exists() and not force:
            raise CertificateSetupError(
                "Refusing to replace the installed MinKNOW client certificate at "
                f"{installed_certificate}. Use --force."
            )

    @staticmethod
    def _run(*command: str) -> subprocess.CompletedProcess[str]:
        try:
            return subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as error:
            detail = error.stderr.strip() or error.stdout.strip()
            raise CertificateSetupError(detail or "OpenSSL command failed") from error
