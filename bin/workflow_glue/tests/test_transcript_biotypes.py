"""Tests for annotation-driven Oarfish transcript-biotype composition."""

import csv
from io import StringIO
from pathlib import Path

import pytest

from workflow_glue import transcript_biotypes as tb
from workflow_glue.oarfish_gtf import convert_lines


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def rows_by_feature(rows):
    return {
        row["feature_id"]: row
        for row in rows
        if row["feature_id"] == row["alias"]
    }


def test_gtf_uses_transcript_fields_then_linked_gene_fallback(tmp_path):
    annotation = write(
        tmp_path / "annotation.gtf",
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; '
        'gene_biotype "protein_coding";\n'
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; '
        'transcript_id "tx1";\n'
        'chr1\ttest\tgene\t200\t250\t.\t+\t.\tgene_id "g2"; '
        'gene_biotype "tRNA";\n'
        'chr1\ttest\ttranscript\t200\t250\t.\t+\t.\tgene_id "g2"; '
        'transcript_id "tx2"; transcript_type "lncRNA";\n'
        'chr1\ttest\ttranscript\t300\t350\t.\t+\t.\tgene_id "g3"; '
        'transcript_id "tx3";\n',
    )

    rows = rows_by_feature(tb.build_biotype_rows(tb.read_annotation(annotation)))

    assert rows["tx1"]["biotype"] == "Protein-coding"
    assert rows["tx2"]["biotype"] == "lncRNA"
    assert rows["tx3"]["biotype"] == "Unknown"


def test_linked_gene_fields_follow_declared_precedence(tmp_path):
    annotation = write(
        tmp_path / "annotation.gtf",
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; '
        'gene_biotype "rRNA"; gene_type "tRNA"; biotype "lncRNA";\n'
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; '
        'transcript_id "tx1";\n',
    )

    rows = rows_by_feature(tb.build_biotype_rows(tb.read_annotation(annotation)))

    assert rows["tx1"]["biotype"] == "rRNA"


def test_gff3_parent_fallback_and_namespace_aliases(tmp_path):
    annotation = write(
        tmp_path / "annotation.gff3",
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:g1;gene_biotype=protein_coding\n"
        "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\t"
        "ID=transcript:tx1;Parent=gene:g1\n"
        "chr1\ttest\tgene\t200\t250\t.\t+\t.\t"
        "ID=gene:g2;gene_biotype=rRNA\n"
        "chr1\ttest\trRNA\t200\t250\t.\t+\t.\t"
        "ID=rna:rrsA;Parent=gene:g2\n",
    )

    rows = tb.build_biotype_rows(tb.read_annotation(annotation))
    by_alias = {row["alias"]: row["biotype"] for row in rows}

    assert by_alias["tx1"] == "Protein-coding"
    assert by_alias["rrsA"] == "rRNA"


def test_converted_prokaryotic_gtf_retains_coding_and_structural_rna(tmp_path):
    source = (
        'chr1\tGenbank\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; '
        'gene_biotype "protein_coding";\n'
        'chr1\tGenbank\tCDS\t1\t97\t.\t+\t0\tgene_id "g1"; '
        'transcript_id "tx1"; protein_id "p1";\n'
        'chr1\tGenbank\tgene\t200\t250\t.\t+\t.\tgene_id "g2"; '
        'gene_biotype "tRNA";\n'
        'chr1\tGenbank\ttranscript\t200\t250\t.\t+\t.\tgene_id "g2"; '
        'transcript_id "tx2"; transcript_biotype "tRNA";\n'
        'chr1\tGenbank\texon\t200\t250\t.\t+\t.\tgene_id "g2"; '
        'transcript_id "tx2";\n'
    )
    converted = StringIO()
    convert_lines(source.splitlines(keepends=True), converted)
    annotation = write(tmp_path / "converted.gtf", converted.getvalue())

    rows = rows_by_feature(tb.build_biotype_rows(tb.read_annotation(annotation)))

    assert rows["tx1"]["biotype"] == "Protein-coding"
    assert rows["tx2"]["biotype"] == "tRNA"


def test_synthetic_eukaryotic_annotation_spans_broad_classes(tmp_path):
    annotation = write(
        tmp_path / "eukaryote.gtf",
        "".join(
            f'chr1\ttest\ttranscript\t{index}\t{index + 20}\t.\t+\t.\t'
            f'transcript_id "{feature_id}"; transcript_biotype "{biotype}";\n'
            for index, (feature_id, biotype) in enumerate(
                [
                    ("coding", "protein_coding"),
                    ("long", "lncRNA"),
                    ("structural", "snoRNA"),
                    ("pseudo", "processed_pseudogene"),
                ],
                start=1,
            )
        ),
    )

    rows = rows_by_feature(tb.build_biotype_rows(tb.read_annotation(annotation)))

    assert rows["coding"]["biotype"] == "Protein-coding"
    assert rows["long"]["biotype"] == "lncRNA"
    assert rows["structural"]["biotype"] == "Other ncRNA"
    assert rows["pseudo"]["biotype"] == "Pseudogene"


def test_gff3_url_decoding_preserves_encoded_commas():
    attributes = tb.parse_attributes("ID=transcript%3Atx1;Note=alpha%2Cbeta")

    assert attributes["ID"] == ("transcript:tx1",)
    assert attributes["Note"] == ("alpha,beta",)


def test_ambiguous_namespace_stripped_alias_becomes_unknown(tmp_path):
    mapping = write(
        tmp_path / "map.tsv",
        "feature_id\talias\traw_biotype\tbiotype\n"
        "transcript:tx1\ttx1\tprotein_coding\tProtein-coding\n"
        "rna:tx1\ttx1\tprotein_coding\tProtein-coding\n",
    )

    exact, aliases = tb.read_biotype_map(mapping)

    assert exact["transcript:tx1"] == "Protein-coding"
    assert aliases["tx1"] == "Unknown"


def test_conflicting_high_priority_biotypes_become_unknown(tmp_path):
    annotation = write(
        tmp_path / "conflict.gtf",
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "tx1"; '
        'transcript_biotype "rRNA";\n'
        'chr1\ttest\texon\t1\t100\t.\t+\t.\ttranscript_id "tx1"; '
        'transcript_biotype "tRNA";\n',
    )

    rows = rows_by_feature(tb.build_biotype_rows(tb.read_annotation(annotation)))

    assert rows["tx1"]["biotype"] == "Unknown"


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("protein_coding", "Protein-coding"),
        ("Mt_rRNA", "rRNA"),
        ("tRNA", "tRNA"),
        ("antisense", "lncRNA"),
        ("miRNA", "Other ncRNA"),
        ("RNase_P_RNA", "Other ncRNA"),
        ("processed_pseudogene", "Pseudogene"),
        ("TEC", "Other"),
        (None, "Unknown"),
    ],
)
def test_canonical_biotype_groups(raw, expected):
    assert tb.canonical_biotype(raw) == expected


def test_summary_uses_estimated_reads_and_includes_unknown(tmp_path):
    annotation = write(
        tmp_path / "annotation.gtf",
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "tx1"; '
        'transcript_biotype "protein_coding";\n'
        'chr1\ttest\ttranscript\t200\t250\t.\t+\t.\ttranscript_id "tx2"; '
        'transcript_biotype "rRNA";\n',
    )
    mapping = tmp_path / "mapping.tsv"
    tb.write_biotype_map(annotation, mapping)
    counts_dir = tmp_path / "counts"
    counts_dir.mkdir()
    write(
        counts_dir / "sample.quant",
        "tname\tlen\tnum_reads\n"
        "tx1\t100\t1.5\n"
        "tx1\t100\t0.5\n"
        "tx2\t50\t1\n"
        "missing\t75\t1\n",
    )
    manifest = write(
        tmp_path / "manifest.tsv",
        "name\tgroup\tcount_file\nS1\tcontrol\tsample.quant\n",
    )

    rows = tb.summarize_biotypes(manifest, counts_dir, mapping)
    by_biotype = {row["biotype"]: row for row in rows}

    assert [row["biotype"] for row in rows] == list(tb.BIOTYPE_ORDER)
    assert float(by_biotype["Protein-coding"]["num_reads"]) == 2
    assert float(by_biotype["Protein-coding"]["fraction"]) == 0.5
    assert float(by_biotype["rRNA"]["fraction"]) == 0.25
    assert float(by_biotype["Unknown"]["fraction"]) == 0.25


def test_zero_total_sample_emits_all_zero_fractions(tmp_path):
    annotation = write(
        tmp_path / "annotation.gtf",
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "tx1"; '
        'transcript_biotype "protein_coding";\n',
    )
    mapping = tmp_path / "mapping.tsv"
    tb.write_biotype_map(annotation, mapping)
    counts_dir = tmp_path / "counts"
    counts_dir.mkdir()
    write(counts_dir / "zero.quant", "tname\tlen\tnum_reads\ntx1\t100\t0\n")
    manifest = write(
        tmp_path / "manifest.tsv",
        "name\tgroup\tcount_file\nzero\tcontrol\tzero.quant\n",
    )

    rows = tb.summarize_biotypes(manifest, counts_dir, mapping)

    assert len(rows) == len(tb.BIOTYPE_ORDER)
    assert all(float(row["fraction"]) == 0 for row in rows)


def test_negative_oarfish_counts_are_rejected(tmp_path):
    counts = write(
        tmp_path / "bad.quant",
        "tname\tlen\tnum_reads\ntx1\t100\t-0.1\n",
    )

    with pytest.raises(ValueError, match="negative num_reads"):
        tb._read_oarfish_counts(counts)


@pytest.mark.parametrize("quantity", ["NaN", "inf", "-inf"])
def test_non_finite_oarfish_counts_are_rejected(tmp_path, quantity):
    counts = write(
        tmp_path / "bad.quant",
        f"tname\tlen\tnum_reads\ntx1\t100\t{quantity}\n",
    )

    with pytest.raises(ValueError, match="non-finite num_reads"):
        tb._read_oarfish_counts(counts)


def test_malformed_oarfish_row_is_rejected(tmp_path):
    counts = write(tmp_path / "bad.quant", "tname\tnum_reads\ntx1\n")

    with pytest.raises(ValueError, match="non-numeric num_reads"):
        tb._read_oarfish_counts(counts)


def test_written_summary_contract(tmp_path):
    annotation = write(
        tmp_path / "annotation.gtf",
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id "tx1"; '
        'transcript_biotype "protein_coding";\n',
    )
    mapping = tmp_path / "mapping.tsv"
    tb.write_biotype_map(annotation, mapping)
    counts_dir = tmp_path / "counts"
    counts_dir.mkdir()
    write(counts_dir / "sample.quant", "tname\tlen\tnum_reads\ntx1\t100\t2\n")
    manifest = write(
        tmp_path / "manifest.tsv",
        "name\tgroup\tcount_file\nS1\tcontrol\tsample.quant\n",
    )
    output = tmp_path / "summary.tsv"

    tb.write_biotype_summary(manifest, counts_dir, mapping, output)

    with output.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == [
            "name",
            "group",
            "biotype",
            "num_reads",
            "fraction",
        ]
        assert len(list(reader)) == len(tb.BIOTYPE_ORDER)
