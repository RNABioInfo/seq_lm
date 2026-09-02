import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from seq_run_manager import main
from seq_run_manager.managers.run_manager import RunManager
from seq_run_manager.managers.sequencing_protocol_manager import SequencingProtocolManager
from seq_run_manager.models.manager_error import ManagerError
from seq_run_manager.models.sample import Sample
from seq_run_manager.models.start_run_config import StartRunConfig
from seq_run_manager.models.stop_acquisition_config import StopAcquisitionConfig


class RunManagerCommandsTest(unittest.TestCase):
    def test_start_maps_alias_and_bam_directory_to_protocol(self):
        connection = MagicMock()
        connection.device.get_flow_cell_info.return_value = SimpleNamespace(
            user_specified_product_code="", product_code="FLO-MIN114"
        )
        connection_manager = MagicMock()
        connection_manager.connect_to_positions.return_value = [connection]
        run_manager = RunManager(connection_manager)
        sample = Sample("sample-a", "control", Path("/output/sample-a"))
        config = SimpleNamespace(samples=[sample], kit="SQK-RNA004")
        protocol = SimpleNamespace(identifier="protocol")

        with (
            patch.object(
                SequencingProtocolManager,
                "get_sequencing_protocol",
                return_value=protocol,
            ),
            patch.object(
                run_manager,
                "_RunManager__resolve_basecalling_settings",
                return_value=("model", 9.0),
            ),
            patch.object(
                SequencingProtocolManager,
                "start_sequencing_protocol",
                return_value="run-id",
            ) as start_protocol,
        ):
            acquisitions = run_manager._RunManager__start_acquisitions(config)

        self.assertEqual(acquisitions[0].sample.id, "sample-a")
        self.assertEqual(start_protocol.call_args.args[2], "sample-a")
        self.assertEqual(start_protocol.call_args.args[3], sample.bam_dir)

    def test_stop_matches_active_run(self):
        connection = MagicMock()
        connection_manager = MagicMock()
        connection_manager.connect_to_all_positions.return_value = [connection]
        run_manager = RunManager(connection_manager)

        with (
            patch.object(
                SequencingProtocolManager,
                "get_currently_active_protocol",
                return_value=SimpleNamespace(run_id="requested"),
            ),
            patch.object(
                SequencingProtocolManager,
                "stop_sequencing_protocol",
                return_value=SimpleNamespace(),
            ) as stop_protocol,
        ):
            run_manager.stop_run("requested")

        stop_protocol.assert_called_once_with(connection, "requested")

    def test_stop_reports_missing_and_failed_runs(self):
        connection_manager = MagicMock()
        connection_manager.connect_to_all_positions.return_value = [MagicMock()]
        run_manager = RunManager(connection_manager)

        with patch.object(
            SequencingProtocolManager,
            "get_currently_active_protocol",
            return_value=SimpleNamespace(run_id="other"),
        ):
            with self.assertRaisesRegex(ManagerError, "specified run ID"):
                run_manager.stop_run("requested")

        with (
            patch.object(
                SequencingProtocolManager,
                "get_currently_active_protocol",
                return_value=SimpleNamespace(run_id="requested"),
            ),
            patch.object(
                SequencingProtocolManager,
                "stop_sequencing_protocol",
                return_value=None,
            ),
        ):
            with self.assertRaisesRegex(ManagerError, "Could not stop protocol"):
                run_manager.stop_run("requested")


class EntrypointLifecycleTest(unittest.TestCase):
    def _connection_values(self, root: Path) -> dict[str, object]:
        paths = [root / name for name in ("client.pem", "key.pem", "ca.crt")]
        for path in paths:
            path.touch()
        return {
            "host": "localhost",
            "port": 9501,
            "client_certificate_path": paths[0],
            "client_private_key_path": paths[1],
            "ca_certificate_path": paths[2],
        }

    def test_stop_disconnects_when_manager_fails(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            config = StopAcquisitionConfig(
                **self._connection_values(Path(temporary_directory)), run_id="run-id"
            )
            connection_manager = MagicMock()
            run_manager = MagicMock()
            run_manager.stop_run.side_effect = ManagerError("failed")

            with (
                patch("seq_run_manager.ArgumentParser.parse_cli_arguments", return_value=config),
                patch("seq_run_manager.check_server"),
                patch(
                    "seq_run_manager.managers.connection_manager.ConnectionManager.connected_with",
                    return_value=connection_manager,
                ),
                patch(
                    "seq_run_manager.managers.run_manager.RunManager",
                    return_value=run_manager,
                ),
            ):
                with self.assertRaises(SystemExit):
                    main([])

            connection_manager.disconnect.assert_called_once_with()

    def test_simulated_start_cleans_up_before_and_after(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = StartRunConfig(
                **self._connection_values(root),
                experiment_id="experiment",
                flow_cell_ids=None,
                position_ids=None,
                kit="kit",
                reference_genome_path=None,
                sampling_regions_path=None,
                adaptive_sampling_mode=None,
                basecall_model=None,
                min_qscore=None,
                output_chunk_size=4000,
                samples=[Sample("sample", "control", root)],
                simulate_run=True,
            )
            connection_manager = MagicMock()
            run_manager = MagicMock()

            with (
                patch("seq_run_manager.ArgumentParser.parse_cli_arguments", return_value=config),
                patch("seq_run_manager.check_server"),
                patch(
                    "seq_run_manager.managers.connection_manager.ConnectionManager.connected_with",
                    return_value=connection_manager,
                ),
                patch(
                    "seq_run_manager.managers.run_manager.RunManager",
                    return_value=run_manager,
                ),
            ):
                main([])

            self.assertEqual(connection_manager.remove_all_simulated_positions.call_count, 2)
            connection_manager.disconnect.assert_called_once_with()

    def test_simulation_cleanup_failure_still_disconnects(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = StartRunConfig(
                **self._connection_values(root),
                experiment_id="experiment",
                flow_cell_ids=None,
                position_ids=None,
                kit="kit",
                reference_genome_path=None,
                sampling_regions_path=None,
                adaptive_sampling_mode=None,
                basecall_model=None,
                min_qscore=None,
                output_chunk_size=4000,
                samples=[Sample("sample", "control", root)],
                simulate_run=True,
            )
            connection_manager = MagicMock()
            connection_manager.remove_all_simulated_positions.side_effect = [
                None,
                ManagerError("cleanup failed"),
            ]

            with (
                patch("seq_run_manager.ArgumentParser.parse_cli_arguments", return_value=config),
                patch("seq_run_manager.check_server"),
                patch(
                    "seq_run_manager.managers.connection_manager.ConnectionManager.connected_with",
                    return_value=connection_manager,
                ),
                patch("seq_run_manager.managers.run_manager.RunManager"),
            ):
                with self.assertRaisesRegex(ManagerError, "cleanup failed"):
                    main([])

            connection_manager.disconnect.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
