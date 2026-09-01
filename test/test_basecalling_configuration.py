import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from minknow_api.tools import protocols

from seq_run_manager.managers.run_manager import RunManager
from seq_run_manager.managers.sequencing_protocol_manager import (
    SequencingProtocolManager,
)
from seq_run_manager.models.run_config import RunConfig
from seq_run_manager.models.sample import Sample


class BasecallingConfigurationTest(unittest.TestCase):
    def _run_config(self, root: Path, **overrides) -> RunConfig:
        credential_paths = [root / name for name in ("client.pem", "key.pem", "ca.crt")]
        for credential_path in credential_paths:
            credential_path.touch()

        values = {
            "host": "localhost",
            "port": 9501,
            "client_certificate_path": str(credential_paths[0]),
            "client_private_key_path": str(credential_paths[1]),
            "ca_certificate_path": str(credential_paths[2]),
            "flow_cell_ids": None,
            "position_ids": None,
            "experiment_id": "experiment",
            "run_id": "1",
            "kit": "SQK-RNA004",
            "reference_genome_path": None,
            "sampling_regions_path": None,
            "adaptive_sampling_mode": None,
            "basecall_model": None,
            "min_qscore": None,
            "output_chunk_size": 4000,
            "samples": [Sample(1, 1, root / "run" / "sample")],
            "simulate_run": False,
        }
        values.update(overrides)
        return RunConfig(**values)

    def test_validation_accepts_model_discovery_and_rejects_legacy_config(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            self._run_config(root).validate()

            with self.assertRaisesRegex(Exception, "Dorado simplex model"):
                self._run_config(
                    root, basecall_model="rna_rp4_130bps_hac_prom.cfg"
                ).validate()

    def test_resolves_default_hac_model_and_its_qscore(self):
        simplex_model = SimpleNamespace(
            name="rna004_130bps_hac@v5.2.0",
            variant="HAC",
            default_q_score_cutoff=9.0,
        )
        configuration = SimpleNamespace(
            kits=["SQK-RNA004"],
            flowcells=["FLO-MIN114"],
            sampling_rate=5000,
            simplex_models=[simplex_model],
        )
        minknow_manager = MagicMock()
        minknow_manager.find_basecall_configurations.return_value = [configuration]
        run_manager = RunManager(SimpleNamespace(manager=minknow_manager))
        position_connection = MagicMock()
        position_connection.device.get_flow_cell_info.return_value = SimpleNamespace(
            user_specified_product_code="", product_code="FLO-MIN114"
        )
        protocol = SimpleNamespace(
            tags={"sample rate": SimpleNamespace(int_value=5000)}
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            config = self._run_config(Path(temporary_directory))
            model_name, min_qscore = run_manager._RunManager__resolve_basecalling_settings(  # type: ignore[attr-defined]
                config, position_connection, protocol, "FLO-MIN114"
            )

        self.assertEqual(model_name, simplex_model.name)
        self.assertEqual(min_qscore, 9.0)
        minknow_manager.find_basecall_configurations.assert_called_once_with(
            flow_cell_product_code="FLO-MIN114",
            sequencing_kit="SQK-RNA004",
            sampling_rate=5000,
            include_outdated=False,
        )

    def test_builds_current_protocol_argument_shapes(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            reference = root / "reference.fasta"
            regions = root / "regions.bed"
            config = self._run_config(
                root,
                reference_genome_path=str(reference),
                sampling_regions_path=str(regions),
                adaptive_sampling_mode="enrich",
            )
            device_connection = MagicMock()
            device_connection.protocol.start_protocol.return_value = "run-id"
            protocol = SimpleNamespace(identifier="protocol-id")

            with patch.object(
                protocols, "make_protocol_arguments", return_value=["--arguments"]
            ) as make_arguments:
                run_id = SequencingProtocolManager.start_sequencing_protocol(
                    device_connection,
                    protocol,
                    "sample-id",
                    root / "run" / "sample",
                    config,
                    "rna004_130bps_hac@v5.2.0",
                    9.0,
                )

        self.assertEqual(run_id, "run-id")
        arguments = make_arguments.call_args.kwargs
        self.assertEqual(
            arguments["basecalling"].simplex_model, "rna004_130bps_hac@v5.2.0"
        )
        self.assertEqual(arguments["basecalling"].min_qscore, 9.0)
        self.assertEqual(arguments["read_until"].reference_files, [str(reference)])
        self.assertEqual(arguments["fastq_arguments"].reads_per_file, 4000)
        self.assertIsNone(arguments["fastq_arguments"].batch_duration)


if __name__ == "__main__":
    unittest.main()
