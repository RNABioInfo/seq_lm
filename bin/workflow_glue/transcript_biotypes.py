"""Map annotation biotypes and summarize Oarfish transcript abundance."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .util import wf_parser


BIOTYPE_ORDER = (
    "Protein-coding",
    "rRNA",
    "tRNA",
    "lncRNA",
    "Other ncRNA",
    "Pseudogene",
    "Other",
    "Unknown",
)
TRANSCRIPT_BIOTYPE_FIELDS = ("transcript_biotype", "transcript_type", "biotype")
GENE_BIOTYPE_FIELDS = ("gene_biotype", "gene_type", "biotype")
TRANSCRIPT_FEATURES = {
    "transcript",
    "mrna",
    "ncrna",
    "lncrna",
    "rrna",
    "trna",
    "tmrna",
    "mirna",
    "snrna",
    "snorna",
}
GTF_ATTRIBUTE_RE = re.compile(
    r'(?P<key>[A-Za-z][A-Za-z0-9_.-]*)\s+"(?P<value>(?:\\.|[^"\\])*)"\s*;?'
)
NAMESPACE_RE = re.compile(
    r"^(?:gene|transcript|rna|mrna|protein|cds|feature):",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class AnnotationRecord:
    """The fields needed to connect an annotation record to a transcript."""

    feature: str
    attributes: dict[str, tuple[str, ...]]
    line_number: int


def _split_gff_values(value: str) -> tuple[str, ...]:
    return tuple(
        unquote(part).strip() for part in value.split(",") if unquote(part).strip()
    )


def parse_attributes(text: str) -> dict[str, tuple[str, ...]]:
    """Parse either quoted GTF attributes or key=value GFF3 attributes."""
    attributes: dict[str, list[str]] = defaultdict(list)
    for match in GTF_ATTRIBUTE_RE.finditer(text):
        value = match.group("value").replace(r'\"', '"').replace(r"\\", "\\")
        if value:
            attributes[match.group("key")].append(value)

    for field in text.split(";"):
        if "=" not in field:
            continue
        key, value = field.strip().split("=", 1)
        if not key:
            continue
        attributes[key].extend(_split_gff_values(value))

    return {key: tuple(dict.fromkeys(values)) for key, values in attributes.items()}


def read_annotation(path: Path) -> list[AnnotationRecord]:
    """Read validated nine-column GTF/GFF3 records."""
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Annotation line {line_number} has {len(fields)} columns; "
                    "expected exactly 9."
                )
            records.append(
                AnnotationRecord(
                    feature=fields[2].lower(),
                    attributes=parse_attributes(fields[8]),
                    line_number=line_number,
                )
            )
    if not records:
        raise ValueError("Annotation contains no data records.")
    return records


def identifier_aliases(identifier: str) -> tuple[str, ...]:
    """Return exact, URL-decoded, and namespace-stripped identifier aliases."""
    decoded = unquote(identifier).strip()
    aliases = [identifier.strip(), decoded]
    aliases.extend(NAMESPACE_RE.sub("", value) for value in tuple(aliases))
    return tuple(dict.fromkeys(value for value in aliases if value))


def _attribute_values(
    attributes: dict[str, tuple[str, ...]],
    fields: tuple[str, ...],
) -> set[str]:
    return {
        value.strip()
        for field in fields
        for value in attributes.get(field, ())
        if value.strip()
    }


def canonical_biotype(raw_biotype: str | None) -> str:
    """Collapse provider-specific annotation labels into stable broad classes."""
    if raw_biotype is None or not raw_biotype.strip():
        return "Unknown"
    normalized = re.sub(r"[^a-z0-9]+", "_", raw_biotype.lower()).strip("_")

    if "pseudogene" in normalized:
        return "Pseudogene"
    if normalized in {"protein_coding", "protein_coding_gene", "coding", "mrna"}:
        return "Protein-coding"
    if normalized in {"rrna", "ribosomal_rna", "mt_rrna", "mitochondrial_rrna"}:
        return "rRNA"
    if normalized in {"trna", "transfer_rna", "mt_trna", "mitochondrial_trna"}:
        return "tRNA"
    if normalized in {
        "lncrna",
        "lincrna",
        "long_non_coding_rna",
        "long_noncoding_rna",
        "antisense",
        "processed_transcript",
        "sense_intronic",
        "sense_overlapping",
        "3prime_overlapping_ncrna",
        "bidirectional_promoter_lncrna",
        "macro_lncrna",
    }:
        return "lncRNA"
    other_ncrna_tokens = (
        "mirna",
        "snrna",
        "snorna",
        "scrna",
        "srp_rna",
        "rnase_p_rna",
        "rnase_mrp_rna",
        "tmrna",
        "misc_rna",
        "small_rna",
        "srna",
        "vault_rna",
        "y_rna",
        "ribozyme",
        "guide_rna",
        "telomerase_rna",
        "ncrna",
        "non_coding_rna",
        "noncoding_rna",
    )
    if any(token in normalized for token in other_ncrna_tokens):
        return "Other ncRNA"
    return "Other"


def _single_or_unknown(values: set[str]) -> tuple[str, str]:
    if len(values) != 1:
        return "", "Unknown"
    raw_biotype = next(iter(values))
    return raw_biotype, canonical_biotype(raw_biotype)


def build_biotype_rows(records: list[AnnotationRecord]) -> list[dict[str, str]]:
    """Build exact transcript IDs and aliases with conservative biotype calls."""
    gene_biotypes: dict[str, dict[int, list[set[str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for record in records:
        if record.feature not in {"gene", "pseudogene"}:
            continue
        identifiers = set(record.attributes.get("gene_id", ())) | set(
            record.attributes.get("ID", ())
        )
        for priority, field in enumerate(GENE_BIOTYPE_FIELDS):
            values = _attribute_values(record.attributes, (field,))
            if values:
                for identifier in identifiers:
                    for alias in identifier_aliases(identifier):
                        gene_biotypes[alias][priority].append(values)

    transcript_candidates: dict[str, dict[int, list[set[str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for record in records:
        transcript_ids = set(record.attributes.get("transcript_id", ()))
        if record.feature in TRANSCRIPT_FEATURES:
            transcript_ids.update(record.attributes.get("ID", ()))
        transcript_ids.discard("")
        if not transcript_ids:
            continue
        for transcript_id in transcript_ids:
            transcript_candidates[transcript_id]

        for priority, field in enumerate(TRANSCRIPT_BIOTYPE_FIELDS):
            values = _attribute_values(record.attributes, (field,))
            if values:
                for transcript_id in transcript_ids:
                    transcript_candidates[transcript_id][priority].append(values)

        offset = len(TRANSCRIPT_BIOTYPE_FIELDS)
        for field_index, field in enumerate(GENE_BIOTYPE_FIELDS):
            values = _attribute_values(record.attributes, (field,))
            if values:
                for transcript_id in transcript_ids:
                    transcript_candidates[transcript_id][offset + field_index].append(
                        values
                    )

        parent_ids = set(record.attributes.get("gene_id", ())) | set(
            record.attributes.get("Parent", ())
        )
        linked_values: dict[int, list[set[str]]] = defaultdict(list)
        for parent_id in parent_ids:
            for alias in identifier_aliases(parent_id):
                for priority, values_at_priority in gene_biotypes.get(
                    alias, {}
                ).items():
                    linked_values[priority].extend(values_at_priority)
        if linked_values:
            for priority, values_at_priority in linked_values.items():
                linked_priority = offset + priority
                for transcript_id in transcript_ids:
                    transcript_candidates[transcript_id][linked_priority].extend(
                        values_at_priority
                    )

    rows = []
    for transcript_id in sorted(transcript_candidates):
        priorities = transcript_candidates[transcript_id]
        selected_values: set[str] = set()
        for priority in sorted(priorities):
            selected_values = set().union(*priorities[priority])
            if selected_values:
                break
        raw_biotype, biotype = _single_or_unknown(selected_values)
        for alias in identifier_aliases(transcript_id):
            rows.append(
                {
                    "feature_id": transcript_id,
                    "alias": alias,
                    "raw_biotype": raw_biotype,
                    "biotype": biotype,
                }
            )
    if not rows:
        raise ValueError("Annotation contains no transcript identifiers.")
    return rows


def write_biotype_map(annotation_path: Path, output_path: Path) -> None:
    rows = build_biotype_rows(read_annotation(annotation_path))
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("feature_id", "alias", "raw_biotype", "biotype"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_biotype_map(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Load exact and alias maps, converting conflicting keys to Unknown."""
    exact_values: dict[str, set[str]] = defaultdict(set)
    alias_values: dict[str, set[str]] = defaultdict(set)
    alias_features: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"feature_id", "alias", "biotype"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Biotype map is missing columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            biotype = row["biotype"]
            if biotype not in BIOTYPE_ORDER:
                raise ValueError(f"Biotype map contains unsupported class {biotype!r}.")
            exact_values[row["feature_id"]].add(biotype)
            alias_values[row["alias"]].add(biotype)
            alias_features[row["alias"]].add(row["feature_id"])

    def collapse(values_by_id: dict[str, set[str]]) -> dict[str, str]:
        return {
            identifier: next(iter(values)) if len(values) == 1 else "Unknown"
            for identifier, values in values_by_id.items()
        }

    alias_map = collapse(alias_values)
    for alias, feature_ids in alias_features.items():
        if len(feature_ids) != 1:
            alias_map[alias] = "Unknown"
    return collapse(exact_values), alias_map


def _resolve_biotype(
    feature_id: str,
    exact_map: dict[str, str],
    alias_map: dict[str, str],
) -> str:
    if feature_id in exact_map:
        return exact_map[feature_id]
    candidates = {
        alias_map[alias]
        for alias in identifier_aliases(feature_id)
        if alias in alias_map
    }
    return next(iter(candidates)) if len(candidates) == 1 else "Unknown"


def _read_oarfish_counts(path: Path) -> dict[str, float]:
    counts: dict[str, float] = defaultdict(float)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"tname", "num_reads"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Oarfish file {path} is missing columns: "
                + ", ".join(sorted(missing))
            )
        for row_number, row in enumerate(reader, start=2):
            feature_id = row["tname"].strip()
            if not feature_id:
                raise ValueError(
                    f"Oarfish file {path} has an empty tname at row {row_number}."
                )
            try:
                value = float(row["num_reads"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Oarfish file {path} has a non-numeric num_reads value at row "
                    f"{row_number}."
                ) from exc
            if not math.isfinite(value):
                raise ValueError(
                    f"Oarfish file {path} has a non-finite num_reads value at row "
                    f"{row_number}."
                )
            if value < 0:
                raise ValueError(
                    f"Oarfish file {path} has a negative num_reads value at row "
                    f"{row_number}."
                )
            counts[feature_id] += value
            if not math.isfinite(counts[feature_id]):
                raise ValueError(
                    f"Oarfish file {path} has an overflowing num_reads total for "
                    f"{feature_id!r}."
                )
    return counts


def summarize_biotypes(
    manifest_path: Path,
    counts_dir: Path,
    mapping_path: Path,
) -> list[dict[str, str]]:
    """Return one fixed-order composition row per sample and canonical class."""
    exact_map, alias_map = read_biotype_map(mapping_path)
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "group", "count_file"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Quantification manifest is missing columns: "
                + ", ".join(sorted(missing))
            )
        manifest_rows = list(reader)
    if not manifest_rows:
        raise ValueError("Quantification manifest contains no samples.")
    sample_keys: set[tuple[str, str]] = set()
    for row_number, sample in enumerate(manifest_rows, start=2):
        empty_fields = [
            field
            for field in required
            if not (sample.get(field) or "").strip()
        ]
        if empty_fields:
            raise ValueError(
                f"Quantification manifest row {row_number} has empty fields: "
                + ", ".join(sorted(empty_fields))
            )
        sample_key = (sample["group"], sample["name"])
        if sample_key in sample_keys:
            raise ValueError(
                "Quantification manifest contains duplicate sample "
                f"{sample['group']}/{sample['name']}."
            )
        sample_keys.add(sample_key)

    output_rows = []
    for sample in manifest_rows:
        counts = _read_oarfish_counts(counts_dir / sample["count_file"])
        by_biotype = {biotype: 0.0 for biotype in BIOTYPE_ORDER}
        for feature_id, value in counts.items():
            by_biotype[_resolve_biotype(feature_id, exact_map, alias_map)] += value
        total = sum(by_biotype.values())
        for biotype in BIOTYPE_ORDER:
            value = by_biotype[biotype]
            output_rows.append(
                {
                    "name": sample["name"],
                    "group": sample["group"],
                    "biotype": biotype,
                    "num_reads": f"{value:.15g}",
                    "fraction": f"{value / total:.15g}" if total > 0 else "0",
                }
            )
    return output_rows


def write_biotype_summary(
    manifest_path: Path,
    counts_dir: Path,
    mapping_path: Path,
    output_path: Path,
) -> None:
    rows = summarize_biotypes(manifest_path, counts_dir, mapping_path)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("name", "group", "biotype", "num_reads", "fraction"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main(args) -> None:
    """Run annotation-map extraction or Oarfish summarization."""
    if args.mode == "map":
        if args.annotation is None:
            raise ValueError("--annotation is required in map mode.")
        write_biotype_map(Path(args.annotation), Path(args.output))
        return
    required = {
        "--manifest": args.manifest,
        "--counts-dir": args.counts_dir,
        "--mapping": args.mapping,
    }
    missing = [option for option, value in required.items() if value is None]
    if missing:
        raise ValueError(
            "Summarize mode requires " + ", ".join(missing) + "."
        )
    write_biotype_summary(
        Path(args.manifest),
        Path(args.counts_dir),
        Path(args.mapping),
        Path(args.output),
    )


def argparser() -> argparse.ArgumentParser:
    parser = wf_parser("transcript_biotypes")
    parser.add_argument("--mode", choices=("map", "summarize"), required=True)
    parser.add_argument("--annotation")
    parser.add_argument("--manifest")
    parser.add_argument("--counts-dir")
    parser.add_argument("--mapping")
    parser.add_argument("--output", required=True)
    return parser
