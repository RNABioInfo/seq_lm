import tempfile
import unittest
from pathlib import Path

import pandas as pd

from deseq_analysis import _load_count_matrix
from deseq_analysis.plots.de_heatmap_plot import create_de_heatmap_plot, z_score
from deseq_analysis.util.metadata import get_metadata, read_quantification_counts


class LiveQuantificationTest(unittest.TestCase):
    def test_manifest_and_oarfish_counts_build_integer_matrix(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            quant_dir = root / "quant"
            quant_dir.mkdir()

            samples = [
                ("rep_1", "CONTROL", [10.4, 20.6]),
                ("rep_2", "control", [11.2, 19.8]),
                ("rep_1", "time_point_1", [30.7, 40.2]),
                ("rep_2", "time_point_1", [29.9, 41.1]),
            ]
            manifest_rows = []
            for index, (name, group, counts) in enumerate(samples, start=1):
                quant_path = quant_dir / f"sample_{index}.quant"
                quant_path.write_text(
                    "tname\tlen\tnum_reads\n"
                    f"tx_1\t100\t{counts[0]}\n"
                    f"tx_2\t200\t{counts[1]}\n"
                )
                manifest_rows.append(f"{name}\t{group}\tquant/{quant_path.name}")

            manifest_path = root / "quant_manifest.tsv"
            manifest_path.write_text(
                "name\tgroup\tcount_file\n" + "\n".join(manifest_rows) + "\n"
            )

            metadata = get_metadata(quant_manifest=str(manifest_path))
            self.assertEqual(metadata["group"].tolist()[:2], ["control", "control"])
            self.assertEqual(len(set(metadata["sample"])), 4)

            matrix = _load_count_matrix(metadata)
            self.assertEqual(matrix.shape, (4, 2))
            self.assertTrue(all(dtype.kind in "iu" for dtype in matrix.dtypes))
            self.assertEqual(matrix.loc["control/rep_1", "tx_1"], 10)
            self.assertEqual(matrix.loc["control/rep_1", "tx_2"], 21)

    def test_featurecounts_input_remains_supported(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            counts_path = Path(temporary_directory) / "counts.txt"
            pd.DataFrame(
                {
                    "Geneid": ["gene_1", "gene_2"],
                    "Length": [100, 200],
                    "sample.bam": [12, 34],
                }
            ).to_csv(counts_path, sep="\t", index=False)

            counts = read_quantification_counts(str(counts_path))
            self.assertEqual(counts.to_dict(), {"gene_1": 12, "gene_2": 34})

    def test_duplicate_group_and_name_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            quant_path = root / "sample.quant"
            quant_path.write_text("tname\tlen\tnum_reads\ntx_1\t100\t30\n")
            manifest_path = root / "quant_manifest.tsv"
            manifest_path.write_text(
                "name\tgroup\tcount_file\n"
                "rep_1\tcontrol\tsample.quant\n"
                "rep_1\tcontrol\tsample.quant\n"
            )

            with self.assertRaisesRegex(ValueError, "duplicate samples"):
                get_metadata(quant_manifest=str(manifest_path))

    def test_heatmap_accepts_oarfish_index_and_one_significant_target(self):
        samples = [
            "control/sample_1",
            "control/sample_2",
            "time_point_1/sample_3",
            "time_point_1/sample_4",
        ]
        counts = pd.DataFrame(
            [[100.0, 110.0, 40.0, 45.0]],
            index=pd.Index(["unassigned_transcript_491"], name="tname"),
            columns=samples,
        )
        metadata = pd.DataFrame(
            {
                "group": pd.Categorical(
                    ["control", "control", "time_point_1", "time_point_1"]
                )
            },
            index=samples,
        )

        standardized = z_score(counts)
        self.assertFalse(standardized.isna().any().any())
        self.assertAlmostEqual(standardized.mean(axis=1).iloc[0], 0.0)
        self.assertIsNotNone(
            create_de_heatmap_plot(
                "time_point_1 vs control",
                counts,
                metadata,
            )
        )


if __name__ == "__main__":
    unittest.main()
