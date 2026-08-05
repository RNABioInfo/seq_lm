"""Convert an NCBI prokaryotic GTF into Oarfish transcript models."""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


ATTRIBUTE_RE = re.compile(
    r'(?P<key>[A-Za-z][A-Za-z0-9_.-]*)\s+"(?P<value>(?:\\.|[^"\\])*)"\s*;'
)


@dataclass(frozen=True)
class GtfRecord:
    """One parsed, non-comment GTF record."""

    fields: tuple[str, ...]
    attributes: dict[str, str]
    line_number: int

    @property
    def feature(self) -> str:
        return self.fields[2]

    @property
    def gene_id(self) -> str | None:
        return self.attributes.get("gene_id")


@dataclass(frozen=True)
class ConversionSummary:
    """Counts used to validate a completed conversion."""

    protein_coding_genes: int
    synthesized_transcripts: int
    synthesized_exons: int
    existing_transcripts: int
    existing_exons: int
    total_transcripts: int
    total_exons: int


def parse_attributes(text: str) -> dict[str, str]:
    """Parse quoted GTF attributes without splitting semicolons inside values."""
    return {
        match.group("key"): match.group("value")
        for match in ATTRIBUTE_RE.finditer(text)
    }


def parse_record(line: str, line_number: int) -> GtfRecord | None:
    """Parse a data line, returning ``None`` for comments and blank lines."""
    if not line.strip() or line.startswith("#"):
        return None
    fields = tuple(line.rstrip("\n").split("\t"))
    if len(fields) != 9:
        raise ValueError(
            f"Line {line_number} has {len(fields)} columns; expected exactly 9."
        )
    try:
        start = int(fields[3])
        end = int(fields[4])
    except ValueError as exc:
        raise ValueError(
            f"Line {line_number} has non-integer GTF coordinates."
        ) from exc
    if start < 1 or end < start:
        raise ValueError(
            f"Line {line_number} has invalid coordinates {start}-{end}."
        )
    return GtfRecord(fields, parse_attributes(fields[8]), line_number)


def _quote_attribute(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _format_attributes(
    gene_attributes: dict[str, str],
    cds_attributes: dict[str, str],
    transcript_id: str,
    *,
    exon: bool,
) -> str:
    values: list[tuple[str, str]] = [
        ("gene_id", gene_attributes["gene_id"]),
        ("transcript_id", transcript_id),
    ]
    for key in ("gene", "locus_tag"):
        value = gene_attributes.get(key)
        if value:
            values.append((key, value))
    protein_id = cds_attributes.get("protein_id")
    if protein_id:
        values.append(("protein_id", protein_id))
    values.append(("transcript_biotype", "protein_coding"))
    if exon:
        values.append(("exon_number", "1"))
    return " ".join(
        f'{key} "{_quote_attribute(value)}";' for key, value in values
    )


def _synthetic_record(
    gene: GtfRecord,
    cds: GtfRecord,
    transcript_id: str,
    feature: str,
) -> str:
    fields = list(gene.fields)
    fields[2] = feature
    fields[5] = "."
    fields[7] = "."
    fields[8] = _format_attributes(
        gene.attributes,
        cds.attributes,
        transcript_id,
        exon=feature == "exon",
    )
    return "\t".join(fields) + "\n"


def _load_records(lines: list[str]) -> list[GtfRecord | None]:
    return [parse_record(line, index) for index, line in enumerate(lines, start=1)]


def _collect_cds_transcripts(
    records: Iterable[GtfRecord | None],
) -> dict[str, dict[str, GtfRecord]]:
    by_gene: dict[str, dict[str, GtfRecord]] = {}
    for record in records:
        if record is None or record.feature != "CDS":
            continue
        gene_id = record.gene_id
        transcript_id = record.attributes.get("transcript_id")
        if not gene_id or not transcript_id:
            raise ValueError(
                f"CDS line {record.line_number} lacks gene_id or transcript_id."
            )
        by_gene.setdefault(gene_id, {}).setdefault(transcript_id, record)
    return by_gene


def convert_lines(lines: list[str], output: TextIO) -> ConversionSummary:
    """Convert GTF lines and write a validated Oarfish-compatible GTF."""
    records = _load_records(lines)
    cds_by_gene = _collect_cds_transcripts(records)
    existing_transcript_ids = {
        record.attributes["transcript_id"]
        for record in records
        if record is not None
        and record.feature == "transcript"
        and record.attributes.get("transcript_id")
    }
    existing_transcripts = sum(
        record is not None and record.feature == "transcript" for record in records
    )
    existing_exons = sum(
        record is not None and record.feature == "exon" for record in records
    )

    protein_coding_gene_ids: set[str] = set()
    synthesized_transcripts = 0
    synthesized_exons = 0
    provenance_written = False

    for line, record in zip(lines, records):
        if not provenance_written and record is not None:
            output.write(
                "#!oarfish-gtf-convert protein-coding genes represented as "
                "single-exon transcripts using NCBI gene spans\n"
            )
            provenance_written = True
        output.write(line if line.endswith("\n") else line + "\n")

        if record is None or record.feature != "gene":
            continue
        if record.attributes.get("gene_biotype") != "protein_coding":
            continue
        gene_id = record.gene_id
        if not gene_id:
            raise ValueError(
                f"Protein-coding gene line {record.line_number} lacks gene_id."
            )
        if gene_id in protein_coding_gene_ids:
            raise ValueError(f"Duplicate protein-coding gene_id {gene_id!r}.")
        protein_coding_gene_ids.add(gene_id)

        transcript_records = cds_by_gene.get(gene_id, {})
        if len(transcript_records) != 1:
            observed = ", ".join(sorted(transcript_records)) or "none"
            raise ValueError(
                f"Protein-coding gene {gene_id!r} must map to exactly one CDS "
                f"transcript_id; observed {observed}."
            )
        transcript_id, cds = next(iter(transcript_records.items()))
        if transcript_id in existing_transcript_ids:
            raise ValueError(
                f"Protein-coding transcript_id {transcript_id!r} already has an "
                "explicit transcript record; refusing to duplicate it."
            )
        output.write(_synthetic_record(record, cds, transcript_id, "transcript"))
        output.write(_synthetic_record(record, cds, transcript_id, "exon"))
        synthesized_transcripts += 1
        synthesized_exons += 1

    if not protein_coding_gene_ids:
        raise ValueError("No protein-coding gene records were found.")

    return ConversionSummary(
        protein_coding_genes=len(protein_coding_gene_ids),
        synthesized_transcripts=synthesized_transcripts,
        synthesized_exons=synthesized_exons,
        existing_transcripts=existing_transcripts,
        existing_exons=existing_exons,
        total_transcripts=existing_transcripts + synthesized_transcripts,
        total_exons=existing_exons + synthesized_exons,
    )


def convert_file(
    input_path: Path,
    output_path: Path,
    *,
    force: bool = False,
) -> ConversionSummary:
    """Convert one file atomically, refusing accidental replacement by default."""
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    if input_path == output_path:
        raise ValueError("Input and output paths must differ.")
    if output_path.exists() and not force:
        raise FileExistsError(
            f"Output already exists: {output_path}. Pass --force to replace it."
        )
    if not input_path.is_file():
        raise FileNotFoundError(f"Input GTF does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = input_path.read_text(encoding="utf-8").splitlines(keepends=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            summary = convert_lines(lines, temporary)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, output_path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="NCBI prokaryotic input GTF")
    parser.add_argument("output", type=Path, help="Oarfish-compatible output GTF")
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing output file",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = convert_file(args.input, args.output, force=args.force)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"protein-coding genes: {summary.protein_coding_genes}")
    print(f"synthesized transcripts: {summary.synthesized_transcripts}")
    print(f"preserved transcripts: {summary.existing_transcripts}")
    print(f"total transcripts: {summary.total_transcripts}")
    print(f"total exons: {summary.total_exons}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
