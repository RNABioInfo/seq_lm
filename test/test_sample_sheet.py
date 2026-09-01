import tempfile
import unittest
from pathlib import Path

from seq_run_manager.utils.sample_sheet import SampleSheetError, parse_sample_sheet


class SampleSheetTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _write(self, contents: str) -> Path:
        samplesheet = self.root / "samples.csv"
        samplesheet.write_text(contents)
        return samplesheet

    def test_parses_optional_fields(self):
        bam_dir = self.root / "bams"
        bam_dir.mkdir()
        samples = parse_sample_sheet(
            self._write(
                "alias,group,bam_dir,is_live,order\n"
                f"sample,control,{bam_dir},false,3\n"
            )
        )
        self.assertEqual(samples[0].alias, "sample")
        self.assertFalse(samples[0].is_live)
        self.assertEqual(samples[0].order, 3)

    def test_rejects_empty_sheet_and_missing_fields(self):
        with self.assertRaisesRegex(SampleSheetError, "at least one row"):
            parse_sample_sheet(self._write("alias,group,bam_dir\n"))
        with self.assertRaisesRegex(SampleSheetError, "missing required fields"):
            parse_sample_sheet(self._write("alias,group\nsample,control\n"))

    def test_rejects_invalid_boolean_order_and_path(self):
        bam_dir = self.root / "bams"
        bam_dir.mkdir()
        with self.assertRaisesRegex(SampleSheetError, "Invalid is_live"):
            parse_sample_sheet(
                self._write(f"alias,group,bam_dir,is_live\nsample,control,{bam_dir},yes\n")
            )
        with self.assertRaisesRegex(SampleSheetError, "Invalid order"):
            parse_sample_sheet(
                self._write(f"alias,group,bam_dir,order\nsample,control,{bam_dir},first\n")
            )
        with self.assertRaisesRegex(SampleSheetError, "does not exist"):
            parse_sample_sheet(
                self._write("alias,group,bam_dir\nsample,control,/missing/bam/directory\n")
            )

    def test_rejects_duplicate_identity_and_start_duplicate_alias(self):
        first = self.root / "first"
        second = self.root / "second"
        first.mkdir()
        second.mkdir()
        with self.assertRaisesRegex(SampleSheetError, "unique within each group"):
            parse_sample_sheet(
                self._write(
                    "alias,group,bam_dir\n"
                    f"sample,control,{first}\n"
                    f"sample,control,{second}\n"
                )
            )

        cross_group = self._write(
            "alias,group,bam_dir\n"
            f"sample,control,{first}\n"
            f"sample,treated,{second}\n"
        )
        self.assertEqual(len(parse_sample_sheet(cross_group)), 2)
        with self.assertRaisesRegex(SampleSheetError, "globally unique"):
            parse_sample_sheet(cross_group, require_unique_aliases=True)


if __name__ == "__main__":
    unittest.main()
