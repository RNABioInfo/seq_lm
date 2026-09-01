import argparse
import csv
from collections.abc import Sequence
from pathlib import Path

from ..models.certificate_setup_config import CertificateSetupConfig
from ..models.run_config import RunConfig
from ..models.sample import Sample


class ArgumentParser:
    @staticmethod
    def parse_cli_arguments(arguments: Sequence[str] | None = None) -> RunConfig:
        parser = argparse.ArgumentParser(
            description="Run manager for MinKNOW runs",
            epilog=(
                "Additional command: seq-run-manager setup-certificates --help"
            ),
        )
        parser.add_argument(
            "--host",
            default="host.docker.internal",
            help="Specify which host to connect to. (Default: localhost)",
        )
        parser.add_argument(
            "--port",
            default=9501,
            help="Specify which port to connect to. (Default: 9501)",
            type=int,
        )
        parser.add_argument(
            "--client_certificate_path",
            "--certificate_path",
            dest="client_certificate_path",
            required=True,
            help="Path to the PEM-encoded client certificate chain",
        )
        parser.add_argument(
            "--client_private_key_path",
            "--key_path",
            dest="client_private_key_path",
            required=True,
            help="Path to the PEM-encoded client private key",
        )
        parser.add_argument(
            "--ca_certificate_path",
            required=True,
            help="Path to the PEM-encoded MinKNOW root CA certificate",
        )
        parser.add_argument(
            "-e", "--experiment_id", required=True, help="Experiment ID (required)"
        )
        parser.add_argument(
            "--run_id", required=True, help="Run number (required)", type=str
        )
        identifier_group = parser.add_mutually_exclusive_group()
        identifier_group.add_argument(
            "-f", "--flow_cell_ids", help="Flow cell ID exclusive to position ID"
        )
        identifier_group.add_argument(
            "-p", "--position_ids", help="Position ID exclusive to flow cell ID"
        )
        parser.add_argument("-k", "--kit", required=True, help="Kit name (required)")
        parser.add_argument("-g", "--reference_genome", help="Reference genome (FASTA)")
        parser.add_argument(
            "-i", "--sampling_regions", help="Regions for adaptive sampling (BED)"
        )
        parser.add_argument(
            "-a",
            "--adaptive_sampling_mode",
            choices=["deplete", "enrich"],
            help="Adaptive sampling mode",
        )
        parser.add_argument(
            "-b",
            "--basecall_model",
            help=(
                "Dorado simplex model name. MinKNOW's default HAC model for the "
                "flow cell, kit, and sampling rate is used when omitted"
            ),
        )
        parser.add_argument(
            "--min_qscore",
            type=float,
            help="Minimum basecall Q-score (default: selected model's cutoff)",
        )
        parser.add_argument(
            "-u",
            "--output_chunk_size",
            help="Output chunk size (default: 4000)",
            default=4000,
            type=int,
        )
        parser.add_argument(
            "-m", "--metadata", help="Metadata file (required)", required=True
        )

        parser.add_argument(
            "--simulate_run", action="store_true", help="Simulate run (default: False)"
        )

        args = parser.parse_args(arguments)

        samples = ArgumentParser.parse_tsv_to_samples(args.metadata)

        config = RunConfig(
            host=args.host,
            port=args.port,
            client_certificate_path=args.client_certificate_path,
            client_private_key_path=args.client_private_key_path,
            ca_certificate_path=args.ca_certificate_path,
            run_id=args.run_id,
            experiment_id=args.experiment_id,
            flow_cell_ids=args.flow_cell_ids,
            position_ids=args.position_ids,
            kit=args.kit,
            reference_genome_path=args.reference_genome,
            sampling_regions_path=args.sampling_regions,
            adaptive_sampling_mode=args.adaptive_sampling_mode,
            basecall_model=args.basecall_model,
            min_qscore=args.min_qscore,
            output_chunk_size=args.output_chunk_size,
            samples=samples,
            simulate_run=args.simulate_run,
        )

        config.validate()
        return config

    @staticmethod
    def parse_certificate_setup_arguments(
        arguments: Sequence[str] | None = None,
    ) -> CertificateSetupConfig:
        parser = argparse.ArgumentParser(
            prog="seq-run-manager setup-certificates",
            description=(
                "Generate MinKNOW client credentials and copy the MinKNOW root CA"
            ),
        )
        parser.add_argument(
            "--output-directory",
            type=Path,
            default=Path.home() / ".config" / "seq-run-manager" / "minknow",
            help="Credential output directory",
        )
        parser.add_argument(
            "--ca-certificate-source",
            type=Path,
            help=(
                "Existing MinKNOW ca.crt; auto-detected from MinKNOW and WSL "
                "locations when omitted"
            ),
        )
        parser.add_argument(
            "--minknow-client-certs-directory",
            type=Path,
            help=(
                "MinKNOW conf/rpc-client-certs directory in which to install the "
                "public client certificate"
            ),
        )
        parser.add_argument(
            "--common-name",
            default="seq-run-manager",
            help="Common Name for the generated client certificate",
        )
        parser.add_argument(
            "--valid-days",
            type=int,
            default=3650,
            help="Client certificate validity in days (default: 3650)",
        )
        parser.add_argument(
            "--key-size",
            type=int,
            choices=(2048, 3072, 4096),
            default=4096,
            help="RSA private key size (default: 4096)",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Replace existing generated and installed credentials",
        )
        args = parser.parse_args(arguments)
        if args.valid_days < 1:
            parser.error("--valid-days must be greater than zero")

        return CertificateSetupConfig(
            output_directory=args.output_directory,
            ca_certificate_source=args.ca_certificate_source,
            minknow_client_certs_directory=args.minknow_client_certs_directory,
            common_name=args.common_name,
            valid_days=args.valid_days,
            key_size=args.key_size,
            force=args.force,
        )

    @staticmethod
    def parse_tsv_to_samples(file_path: str) -> list[Sample]:
        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            return [
                Sample(
                    run_number=int(row["run_number"]),
                    replicate_number=int(row["replicate_number"]),
                    replicate_dir=Path(row["replicate_dir"]),
                )
                for row in reader
            ]
