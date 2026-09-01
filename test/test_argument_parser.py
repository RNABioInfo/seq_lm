import tempfile
import unittest
from pathlib import Path

from seq_run_manager.models.certificate_config import CertificateConfig
from seq_run_manager.models.start_run_config import StartRunConfig
from seq_run_manager.models.stop_acquisition_config import StopAcquisitionConfig
from seq_run_manager.utils.argument_parser import ArgumentParser


class ArgumentParserTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.credentials = [self.root / name for name in ("client.pem", "key.pem", "ca.crt")]
        for credential in self.credentials:
            credential.touch()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _connection_arguments(self) -> list[str]:
        return [
            "--client-certificate-path",
            str(self.credentials[0]),
            "--client-private-key-path",
            str(self.credentials[1]),
            "--ca-certificate-path",
            str(self.credentials[2]),
        ]

    def test_parses_certificate_command(self):
        config = ArgumentParser.parse_cli_arguments(["cert", "--valid-days", "30"])
        self.assertIsInstance(config, CertificateConfig)
        self.assertEqual(config.valid_days, 30)

    def test_parses_start_command_and_identifier_list(self):
        first = self.root / "first"
        second = self.root / "second"
        first.mkdir()
        second.mkdir()
        samplesheet = self.root / "samples.csv"
        samplesheet.write_text(
            "alias,group,bam_dir,is_live,order\n"
            f"first,control,{first},true,0\n"
            f"second,treated,{second},false,1\n"
        )

        config = ArgumentParser.parse_cli_arguments(
            [
                "start",
                *self._connection_arguments(),
                "--samplesheet",
                str(samplesheet),
                "--experiment-id",
                "experiment",
                "--kit",
                "SQK-RNA004",
                "--position-ids",
                "P1",
                "P2",
            ]
        )

        self.assertIsInstance(config, StartRunConfig)
        self.assertEqual(config.position_ids, ["P1", "P2"])
        self.assertEqual([sample.id for sample in config.samples], ["first", "second"])
        self.assertEqual(config.samples[0].bam_dir, first)

    def test_parses_stop_command(self):
        config = ArgumentParser.parse_cli_arguments(
            ["stop", *self._connection_arguments(), "--run-id", "run-123"]
        )
        self.assertIsInstance(config, StopAcquisitionConfig)
        self.assertEqual(config.run_id, "run-123")

    def test_requires_subcommand_and_rejects_removed_command(self):
        with self.assertRaises(SystemExit):
            ArgumentParser.parse_cli_arguments([])
        with self.assertRaises(SystemExit):
            ArgumentParser.parse_cli_arguments(["setup-certificates"])

    def test_flow_cell_and_position_ids_are_mutually_exclusive(self):
        with self.assertRaises(SystemExit):
            ArgumentParser.build_parser().parse_args(
                [
                    "start",
                    *self._connection_arguments(),
                    "--samplesheet",
                    "samples.csv",
                    "--experiment-id",
                    "experiment",
                    "--kit",
                    "kit",
                    "--position-ids",
                    "P1",
                    "--flow-cell-ids",
                    "F1",
                ]
            )


if __name__ == "__main__":
    unittest.main()
