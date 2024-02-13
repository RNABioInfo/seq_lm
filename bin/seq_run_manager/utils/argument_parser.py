import argparse
import csv
from pathlib import Path
from typing import List

from ..models.run_config import RunConfig
from ..models.sample import Sample


class ArgumentParser:
    @staticmethod
    def parse_cli_arguments() -> RunConfig:
        parser = argparse.ArgumentParser(description="Run manager for MinKNOW runs")
        parser.add_argument(
            "--host",
            default="localhost",
            help="Specify which host to connect to. (Default: localhost)",
        )
        parser.add_argument(
            "--port",
            default=9501,
            help="Specify which port to connect to. (Default: 9501)",
        )
        parser.add_argument(
            "--certificate_path",
            required=True,
            help="Specify the path to the certificate (required)",
        )
        parser.add_argument(
            "--key_path", required=True, help="Specify the path to the key (required)"
        )
        parser.add_argument(
            "--replicate_count",
            required=True,
            help="Count of replicates in run (required)",
            type=int,
        )
        parser.add_argument(
            "--run_number", required=True, help="Run number (required)", type=int
        )
        identifier_group = parser.add_mutually_exclusive_group()
        identifier_group.add_argument(
            "-f", "--flow_cell_ids", help="Flow cell ID exclusive to position ID"
        )
        identifier_group.add_argument(
            "-p", "--position_ids", help="Position ID exclusive to flow cell ID"
        )
        parser.add_argument(
            "-e", "--experiment_id", required=True, help="Experiment ID (required)"
        )
        parser.add_argument("-k", "--kit", required=True, help="Kit name (required)")
        parser.add_argument(
            "-g", "--reference_genome", required=True, help="Reference genome (FASTA)"
        )
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
            "--basecall_config",
            help="Basecall config (default: rna_rp4_130bps_sup_prom.cfg)",
            default="rna_rp4_130bps_hac_prom",
        )
        parser.add_argument(
            "-u",
            "--output_chunk_size",
            help="Output chunk size (default: 4000)",
            default=4000,
        )
        parser.add_argument(
            "-m", "--metadata", help="Metadata file (required)", required=True
        )

        parser.add_argument(
            "--simulate_run", action="store_true", help="Simulate run (default: False)"
        )

        args = parser.parse_args()

        samples = ArgumentParser.parse_tsv_to_samples(args.metadata)

        return RunConfig(
            host=args.host,
            port=args.port,
            certificate_path=args.certificate_path,
            key_path=args.key_path,
            run_number=args.run_number,
            replicate_count=args.replicate_count,
            flow_cell_ids=args.flow_cell_ids,
            position_ids=args.position_ids,
            experiment_id=args.experiment_id,
            kit=args.kit,
            reference_genome_path=args.reference_genome,
            sampling_regions_path=args.sampling_regions,
            adaptive_sampling_mode=args.adaptive_sampling_mode,
            basecall_config=args.basecall_config,
            output_chunk_size=args.output_chunk_size,
            samples=samples,
            simulate_run=args.simulate_run,
        )

    @staticmethod
    def parse_tsv_to_samples(file_path: str) -> List[Sample]:
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
