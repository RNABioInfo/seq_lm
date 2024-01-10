import argparse
import models.run_config as rc


class ArgumentParser:
    @staticmethod
    def parse_cli_arguments() -> rc.RunConfig:
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
            help="Specify the path to the certificate.",
        )
        parser.add_argument(
            "--key_path", required=True, help="Specify the path to the key."
        )
        identifier_group = parser.add_mutually_exclusive_group(required=True)
        identifier_group.add_argument(
            "-f", "--flow_cell_id", help="Flow cell ID exclusive to position ID"
        )
        identifier_group.add_argument(
            "-p", "--position_id", help="Position ID exclusive to flow cell ID"
        )
        parser.add_argument("-s", "--sample_id", required=True, help="Sample ID")
        parser.add_argument(
            "-e", "--experiment_id", required=True, help="Experiment ID"
        )
        parser.add_argument("-k", "--kit", required=True, help="Kit name")
        parser.add_argument(
            "-r",
            "--flowcell_product_code",
            required=True,
            help="Flow cell product code",
        )
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
            help="Basecall config",
            default="rna_rp4_130bps_sup_prom.cfg",
        )
        parser.add_argument(
            "-u", "--output_chunk_size", help="Output chunk size", default=2000
        )
        parser.add_argument("-o", "--output_dir", help="Output directory", default=".")

        args = parser.parse_args()

        return rc.RunConfig(
            host=args.host,
            port=args.port,
            certificate_path=args.certificate_path,
            key_path=args.key_path,
            flow_cell_id=args.flow_cell_id,
            position_id=args.position_id,
            sample_id=args.sample_id,
            experiment_id=args.experiment_id,
            flow_cell_product_code=args.flowcell_product_code,
            kit=args.kit,
            reference_genome_path=args.reference_genome,
            sampling_regions_path=args.sampling_regions,
            adaptive_sampling_mode=args.adaptive_sampling_mode,
            basecall_config=args.basecall_config,
            output_chunk_size=args.output_chunk_size,
            output_dir=args.output_dir,
        )
