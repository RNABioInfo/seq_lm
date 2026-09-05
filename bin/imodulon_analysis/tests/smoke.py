"""Dependency and CLI smoke test; also shared by integration fixtures."""

import json
from pathlib import Path
import subprocess
import tempfile


def fixture(root):
    root = Path(root)
    (root / "matrix.csv").write_text(",positive,negative\ng1,1,0\ng2,0,-2\n")
    (root / "annotation.gtf").write_text(
        'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; locus_tag "LOC1";\n'
        'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1b";\n'
        'chr\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n'
    )
    rows = ["name\tgroup\torder\tcount_file\tsource_batch_index"]
    for i, values in enumerate(
        [(20, 5, 40, 35), (25, 5, 35, 35), (50, 10, 15, 25), (55, 10, 10, 25)]
    ):
        (root / f"q{i}.quant").write_text(
            "tname\tlen\tnum_reads\n"
            + "".join(
                f"{t}\t100\t{n}\n" for t, n in zip(["t1", "t1b", "t2", "other"], values)
            )
        )
        rows.append(
            f"rep{i % 2}\t{'control' if i < 2 else 'treated'}\t{0 if i < 2 else 10}\tq{i}.quant\t0"
        )
    (root / "manifest.tsv").write_text("\n".join(rows) + "\n")
    return root


def main():
    with tempfile.TemporaryDirectory(prefix="ica smoke ") as tmp:
        root = fixture(tmp)
        subprocess.run(
            [
                "imodulon-analysis",
                "prepare",
                "--matrix",
                str(root / "matrix.csv"),
                "--annotation",
                str(root / "annotation.gtf"),
                "--output",
                str(root / "prepared"),
            ],
            check=True,
        )
        subprocess.run(
            [
                "imodulon-analysis",
                "analyze",
                "--prepared",
                str(root / "prepared"),
                "--manifest",
                str(root / "manifest.tsv"),
                "--counts-dir",
                str(root),
                "--min-reads",
                "0",
                "--output",
                str(root / "results"),
            ],
            check=True,
        )
        assert (
            json.loads((root / "results/status.json").read_text())["tested_count"] == 2
        )
        assert (root / "results/differential_activity.tsv").exists()
    print("iModulon CLI smoke test passed")


if __name__ == "__main__":
    main()
