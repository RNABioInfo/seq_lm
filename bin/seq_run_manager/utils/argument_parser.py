import argparse
from collections.abc import Sequence
from pathlib import Path

from ..models.certificate_config import CertificateConfig
from ..models.manager_error import ManagerError
from ..models.start_run_config import StartRunConfig
from ..models.stop_acquisition_config import StopAcquisitionConfig
from .sample_sheet import SampleSheetError, parse_sample_sheet

CommandConfig = CertificateConfig | StartRunConfig | StopAcquisitionConfig


class ArgumentParser:
    @staticmethod
    def _add_connection_arguments(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--host",
            default="host.docker.internal",
            help="MinKNOW host (default: host.docker.internal)",
        )
        parser.add_argument(
            "--port", default=9501, type=int, help="MinKNOW manager port (default: 9501)"
        )
        parser.add_argument(
            "--client-certificate-path",
            required=True,
            type=Path,
            help="Path to the PEM-encoded client certificate chain",
        )
        parser.add_argument(
            "--client-private-key-path",
            required=True,
            type=Path,
            help="Path to the PEM-encoded client private key",
        )
        parser.add_argument(
            "--ca-certificate-path",
            required=True,
            type=Path,
            help="Path to the PEM-encoded MinKNOW root CA certificate",
        )

    @staticmethod
    def _add_certificate_arguments(parser: argparse.ArgumentParser) -> None:
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

    @staticmethod
    def build_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            prog="seq-run-manager",
            description="Manage MinKNOW certificates and sequencing acquisitions",
        )
        subparsers = parser.add_subparsers(dest="command", required=True)

        certificate_parser = subparsers.add_parser(
            "cert",
            help="Generate MinKNOW client credentials",
            description="Generate MinKNOW client credentials and copy the MinKNOW root CA",
        )
        ArgumentParser._add_certificate_arguments(certificate_parser)

        start_parser = subparsers.add_parser(
            "start", help="Start acquisitions from a workflow samplesheet"
        )
        ArgumentParser._add_connection_arguments(start_parser)
        start_parser.add_argument(
            "--samplesheet", required=True, type=Path, help="Workflow samplesheet CSV"
        )
        start_parser.add_argument(
            "--experiment-id", required=True, help="MinKNOW protocol group ID"
        )
        identifiers = start_parser.add_mutually_exclusive_group()
        identifiers.add_argument(
            "--flow-cell-ids", nargs="+", help="Flow-cell IDs in samplesheet order"
        )
        identifiers.add_argument(
            "--position-ids", nargs="+", help="Position IDs in samplesheet order"
        )
        start_parser.add_argument("--kit", required=True, help="Sequencing kit name")
        start_parser.add_argument(
            "--reference-genome", type=Path, help="Reference genome (FASTA or MMI)"
        )
        start_parser.add_argument(
            "--sampling-regions", type=Path, help="Adaptive-sampling regions (BED)"
        )
        start_parser.add_argument(
            "--adaptive-sampling-mode", choices=("deplete", "enrich")
        )
        start_parser.add_argument(
            "--basecall-model",
            help="Dorado simplex model name (default: MinKNOW HAC model)",
        )
        start_parser.add_argument(
            "--min-qscore",
            type=float,
            help="Minimum basecall Q-score (default: selected model's cutoff)",
        )
        start_parser.add_argument(
            "--output-chunk-size",
            default=4000,
            type=int,
            help="Reads per output file (default: 4000)",
        )
        start_parser.add_argument(
            "--simulate-run", action="store_true", help="Use simulated positions"
        )

        stop_parser = subparsers.add_parser(
            "stop", help="Stop an active acquisition by run ID"
        )
        ArgumentParser._add_connection_arguments(stop_parser)
        stop_parser.add_argument(
            "--run-id", required=True, help="MinKNOW protocol/acquisition run ID"
        )
        return parser

    @staticmethod
    def _connection_values(args: argparse.Namespace) -> dict[str, object]:
        return {
            "host": args.host,
            "port": args.port,
            "client_certificate_path": args.client_certificate_path,
            "client_private_key_path": args.client_private_key_path,
            "ca_certificate_path": args.ca_certificate_path,
        }

    @staticmethod
    def parse_cli_arguments(
        arguments: Sequence[str] | None = None,
    ) -> CommandConfig:
        parser = ArgumentParser.build_parser()
        args = parser.parse_args(arguments)

        if args.command == "cert":
            if args.valid_days < 1:
                parser.error("--valid-days must be greater than zero")
            return CertificateConfig(
                output_directory=args.output_directory,
                ca_certificate_source=args.ca_certificate_source,
                minknow_client_certs_directory=args.minknow_client_certs_directory,
                common_name=args.common_name,
                valid_days=args.valid_days,
                key_size=args.key_size,
                force=args.force,
            )

        connection_values = ArgumentParser._connection_values(args)
        if args.command == "start":
            try:
                samples = parse_sample_sheet(
                    args.samplesheet, require_unique_aliases=True
                )
            except SampleSheetError as error:
                parser.error(str(error))
            config: StartRunConfig | StopAcquisitionConfig = StartRunConfig(
                **connection_values, # type: ignore
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
        else:
            config = StopAcquisitionConfig(**connection_values, run_id=args.run_id) # type: ignore

        try:
            config.validate()
        except ManagerError as error:
            parser.error(str(error))
        return config
