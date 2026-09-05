"""Prepare a fixed ICA basis and analyze cumulative Oarfish snapshots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

from workflow_glue.transcript_biotypes import (
    TRANSCRIPT_FEATURES,
    identifier_aliases,
    read_annotation,
)
from . import __version__


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_table(path):
    """Read exported CSV/TSV without pandas changing duplicate/blank headers."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        delimiter = "\t" if "\t" in first else ","
        handle.seek(0)
        rows = list(csv.reader(handle, delimiter=delimiter))
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(
            f"{path}: expected a header and data, with at least two columns"
        )
    width = len(rows[0])
    if any(len(row) != width for row in rows[1:]):
        raise ValueError(f"{path}: inconsistent column counts")
    headers = [x.strip() for x in rows[0]]
    if len(set(headers)) != len(headers):
        raise ValueError(f"{path}: duplicate column identifiers")
    return pd.DataFrame(rows[1:], columns=headers)


def identifiers(values, label):
    result = [str(x).strip() for x in values]
    if any(not x or any(c in x for c in "\t\r\n") for x in result):
        raise ValueError(
            f"{label}: identifiers must be nonempty and contain no tabs/newlines"
        )
    if len(result) != len(set(result)):
        raise ValueError(f"{label}: duplicate identifiers")
    return result


def numeric(values, label, nonnegative=False):
    try:
        result = np.asarray(values, dtype=float)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{label}: nonnumeric values") from error
    if not np.isfinite(result).all() or (nonnegative and (result < 0).any()):
        raise ValueError(
            f"{label}: values must be finite"
            + (" and nonnegative" if nonnegative else "")
        )
    return result


def annotation_targets(path):
    """Resolve canonical genes and transcript parentage from GTF or GFF3."""
    records = read_annotation(Path(path))
    gene_aliases = defaultdict(set)
    for record in records:
        a = record.attributes
        ids = a.get("gene_id", ())
        if record.feature in {"gene", "pseudogene"}:
            ids = ids or a.get("ID", ())
        if len(ids) > 1:
            raise ValueError(
                f"Ambiguous gene identifier at annotation line {record.line_number}"
            )
        if ids:
            gene = ids[0]
            aliases = set(ids)
            for field in ("locus_tag", "old_locus_tag"):
                aliases.update(a.get(field, ()))
            if record.feature in {"gene", "pseudogene"}:
                aliases.update(a.get("ID", ()))
            for alias in aliases:
                for variant in identifier_aliases(alias):
                    gene_aliases[variant].add(gene)

    transcript_genes = defaultdict(set)
    transcript_aliases = defaultdict(set)
    for record in records:
        a = record.attributes
        transcripts = set(a.get("transcript_id", ()))
        if not transcripts and record.feature in TRANSCRIPT_FEATURES:
            transcripts.update(a.get("ID", ()))
        if not transcripts:
            continue
        parents = set(a.get("gene_id", ()))
        if record.feature in TRANSCRIPT_FEATURES:
            parents.update(a.get("Parent", ()))
        genes = set()
        for parent in parents:
            for alias in identifier_aliases(parent):
                genes.update(gene_aliases.get(alias, ()))
        if len(genes) != 1:
            raise ValueError(
                f"Missing or ambiguous transcript parentage at annotation line {record.line_number}"
            )
        gene = next(iter(genes))
        for field in ("locus_tag", "old_locus_tag"):
            for alias in a.get(field, ()):
                gene_aliases[alias].add(gene)
        for transcript in transcripts:
            transcript_genes[transcript].add(gene)
            aliases = {transcript}
            if record.feature in TRANSCRIPT_FEATURES:
                aliases.update(a.get("ID", ()))
            for alias in aliases:
                for variant in identifier_aliases(alias):
                    transcript_aliases[variant].add(transcript)
    if not transcript_genes or any(len(g) != 1 for g in transcript_genes.values()):
        raise ValueError(
            "Annotation has no transcripts or ambiguous transcript parentage"
        )
    transcripts = {t: next(iter(g)) for t, g in transcript_genes.items()}
    return gene_aliases, transcripts, transcript_aliases


def prepare(matrix, annotation, gene_map, min_coverage, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    table = read_table(matrix)
    genes = identifiers(table.iloc[:, 0], "Matrix genes")
    components = identifiers(table.columns[1:], "Matrix components")
    weights = numeric(table.iloc[:, 1:], "Matrix weights")
    aliases, transcripts, transcript_aliases = annotation_targets(annotation)
    genes_with_targets = set(transcripts.values())
    explicit = {}
    if gene_map:
        mapping = read_table(gene_map)
        if not {"gene_id", "model_gene_id"}.issubset(mapping.columns):
            raise ValueError("Gene map requires gene_id and model_gene_id columns")
        sources = identifiers(mapping.gene_id, "Map gene_id")
        destinations = identifiers(mapping.model_gene_id, "Map model_gene_id")
        for source, destination in zip(sources, destinations):
            if source not in genes_with_targets or destination not in genes:
                raise ValueError(
                    f"Unknown gene in explicit mapping: {source} -> {destination}"
                )
            explicit[destination] = source
    rows, matches = [], {}
    for gene in genes:
        candidates = (
            ({explicit[gene]} if gene in explicit else set())
            if gene_map
            else aliases.get(gene, set())
        )
        candidates = set(candidates) & genes_with_targets
        if len(candidates) > 1:
            rows.append(
                dict(
                    model_gene_id=gene,
                    gene_id="",
                    transcript_id="",
                    method="alias",
                    status="ambiguous",
                )
            )
        elif candidates:
            source = next(iter(candidates))
            if source in matches.values():
                raise ValueError(
                    f"Annotation gene {source} maps to multiple model genes"
                )
            matches[gene] = source
            for transcript in sorted(t for t, g in transcripts.items() if g == source):
                rows.append(
                    dict(
                        model_gene_id=gene,
                        gene_id=source,
                        transcript_id=transcript,
                        method="explicit" if gene_map else "annotation",
                        status="mapped",
                    )
                )
        else:
            rows.append(
                dict(
                    model_gene_id=gene,
                    gene_id="",
                    transcript_id="",
                    method="explicit" if gene_map else "annotation",
                    status="missing",
                )
            )
    pd.DataFrame(rows).to_csv(output / "gene_mapping.tsv", sep="\t", index=False)
    if any(r["status"] == "ambiguous" for r in rows):
        raise ValueError(
            "Ambiguous model gene aliases; supply an explicit gene map (see gene_mapping.tsv)"
        )
    coverage = len(matches) / len(genes)
    if coverage < min_coverage or not matches:
        raise ValueError(
            f"Model gene coverage {len(matches)}/{len(genes)} ({coverage:.3%}) below required {min_coverage:.3%}; inspect gene_mapping.tsv or explicitly lower ica_min_gene_coverage"
        )
    shared = [g for g in genes if g in matches]
    m = weights[[i for i, g in enumerate(genes) if g in matches]]
    u, singular, vt = np.linalg.svd(m, full_matrices=False)
    tolerance = max(m.shape) * np.finfo(float).eps * singular[0]
    if not np.isfinite(singular).all() or not math.isfinite(tolerance):
        raise ValueError("Matrix scale exceeds numerical precision for SVD")
    rank = int((singular > tolerance).sum())
    diagnostics = dict(
        rank=rank,
        singular_values=singular.tolist(),
        rank_tolerance=float(tolerance),
        condition_number=float(singular[0] / singular[-1])
        if singular[-1] > 0
        else None,
        shared_gene_count=len(shared),
        model_gene_count=len(genes),
        gene_coverage=coverage,
    )
    write_json(output / "diagnostics.json", diagnostics)
    if rank != len(components):
        raise ValueError(
            f"Shared ICA matrix is rank deficient: rank {rank}, components {len(components)}"
        )
    inverse = (vt.T / singular) @ u.T
    if not np.isfinite(inverse).all():
        raise ValueError("Matrix scale produces a nonfinite pseudoinverse")
    transcript_map = {
        r["transcript_id"]: r["model_gene_id"] for r in rows if r["status"] == "mapped"
    }
    ambiguous = [
        alias
        for alias, ts in transcript_aliases.items()
        if len(ts) > 1 and any(t in transcript_map for t in ts)
    ]
    if ambiguous:
        raise ValueError(f"Ambiguous transcript aliases: {ambiguous[:5]}")
    # Ratios are invariant to column scale; normalize before squaring to avoid
    # overflow for otherwise valid, arbitrarily scaled ICA components.
    column_scale = np.max(np.abs(weights), axis=0)
    pd.DataFrame(
        dict(
            component_id=components,
            gene_coverage=coverage,
            retained_squared_weight_fraction=np.sum((m / column_scale) ** 2, axis=0)
            / np.sum((weights / column_scale) ** 2, axis=0),
        )
    ).to_csv(output / "component_coverage.tsv", sep="\t", index=False)
    np.savez(output / "basis.npz", weights=m, inverse=inverse)
    write_json(
        output / "model.json",
        dict(
            schema_version=1,
            genes=shared,
            components=components,
            transcript_map=transcript_map,
            transcript_aliases={
                a: next(iter(ts))
                for a, ts in transcript_aliases.items()
                if len(ts) == 1
            },
            diagnostics=diagnostics,
            min_gene_coverage=min_coverage,
            hashes=dict(
                matrix=sha256(matrix),
                annotation=sha256(annotation),
                gene_map=sha256(gene_map) if gene_map else None,
            ),
        ),
    )


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1)
    return result


STAT_COLUMNS = [
    "component_id",
    "target_group",
    "control_group",
    "activity_difference",
    "target_mean",
    "control_mean",
    "target_sd",
    "control_sd",
    "target_n",
    "control_n",
    "standard_error",
    "degrees_of_freedom",
    "t_statistic",
    "ci_lower",
    "ci_upper",
    "p_value",
    "adjusted_p_value",
    "significant",
    "status",
]


def differential_activity(
    uncentered, centered, metadata, components, control_group, cutoff
):
    rows = []
    control = metadata.group.eq(control_group).to_numpy()
    for group in dict.fromkeys(metadata.group):
        if group == control_group:
            continue
        target = metadata.group.eq(group).to_numpy()
        contrast_rows = []
        for i, component in enumerate(components):
            x, y = uncentered[i, target], uncentered[i, control]
            # Identical floating-point values must stay exactly zero-variance;
            # subtracting their rounded mean can otherwise introduce tiny noise.
            vx = (
                (0.0 if np.ptp(x) == 0 else float(x.var(ddof=1)))
                if len(x) > 1
                else None
            )
            vy = (
                (0.0 if np.ptp(y) == 0 else float(y.var(ddof=1)))
                if len(y) > 1
                else None
            )
            row = dict.fromkeys(STAT_COLUMNS, None)
            row.update(
                component_id=component,
                target_group=group,
                control_group=control_group,
                activity_difference=float(x.mean() - y.mean()),
                target_mean=float(centered[i, target].mean()),
                control_mean=float(centered[i, control].mean()),
                target_sd=math.sqrt(vx) if vx is not None else None,
                control_sd=math.sqrt(vy) if vy is not None else None,
                target_n=len(x),
                control_n=len(y),
                status="insufficient_replicates",
            )
            if min(len(x), len(y)) >= 2:
                variance = vx / len(x) + vy / len(y)
                if not math.isfinite(variance):
                    raise ValueError(
                        "Activity variance is nonfinite; check matrix scaling"
                    )
                if variance == 0:
                    row["status"] = "zero_variance"
                else:
                    se = math.sqrt(variance)
                    df = 1 / (
                        (vx / len(x) / variance) ** 2 / (len(x) - 1)
                        + (vy / len(y) / variance) ** 2 / (len(y) - 1)
                    )
                    t = row["activity_difference"] / se
                    margin = float(stats.t.ppf(0.975, df)) * se
                    row.update(
                        standard_error=se,
                        degrees_of_freedom=df,
                        t_statistic=t,
                        ci_lower=row["activity_difference"] - margin,
                        ci_upper=row["activity_difference"] + margin,
                        p_value=float(2 * stats.t.sf(abs(t), df)),
                        status="tested",
                    )
            contrast_rows.append(row)
        adjusted = bh_adjust(
            [r["p_value"] if r["status"] == "tested" else 1 for r in contrast_rows]
        )
        for row, q in zip(contrast_rows, adjusted):
            if row["status"] == "tested":
                row.update(adjusted_p_value=float(q), significant=bool(q <= cutoff))
        rows.extend(contrast_rows)
    return pd.DataFrame(rows, columns=STAT_COLUMNS)


def stable_sample_id(group, alias):
    return (
        "sample_"
        + hashlib.sha256(
            json.dumps([group, alias], ensure_ascii=False).encode()
        ).hexdigest()
    )


def analyze(
    prepared,
    manifest,
    counts_dir,
    output,
    log_base=2.0,
    pseudocount=1.0,
    min_reads=10000,
    cutoff=0.05,
    batch_index=0,
    analysis_index=0,
    report_sequence=0,
):
    prepared, output = Path(prepared), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    model = json.loads((prepared / "model.json").read_text())
    with np.load(prepared / "basis.npz", allow_pickle=False) as basis:
        m, inverse = basis["weights"], basis["inverse"]
    metadata = read_table(manifest)
    if not {"name", "group", "count_file"}.issubset(metadata.columns):
        raise ValueError("Manifest requires name, group, count_file")
    for column in ("name", "group"):
        metadata[column] = metadata[column].str.strip()
        if metadata[column].eq("").any():
            raise ValueError(f"Missing manifest {column}")
    metadata["sample_id"] = [
        stable_sample_id(g, n) for g, n in zip(metadata.group, metadata.name)
    ]
    identifiers(metadata.sample_id, "Samples")
    control_labels = set(
        metadata.loc[metadata.group.str.lower().eq("control"), "group"]
    )
    if len(control_labels) != 1:
        raise ValueError("Require one unambiguous case-insensitive control group label")
    control_group = next(iter(control_labels))
    control = metadata.group.eq(control_group).to_numpy()
    if control.sum() < 2:
        raise ValueError("Require at least two control samples")
    for field in ("order", "source_batch_index"):
        if field not in metadata:
            metadata[field] = ""
    metadata = metadata.rename(columns={"name": "alias"})
    abundance = np.zeros((len(model["genes"]), len(metadata)))
    positions = {g: i for i, g in enumerate(model["genes"])}
    hashes, totals = [], []
    for sample_index, row in metadata.iterrows():
        path = Path(counts_dir) / row.count_file
        counts = read_table(path)
        if not {"tname", "num_reads"}.issubset(counts.columns):
            raise ValueError(f"{path}: requires tname and num_reads")
        ids = identifiers(counts.tname, f"{path} targets")
        values = numeric(counts.num_reads, f"{path} abundance", nonnegative=True)
        total = float(values.sum())
        if not math.isfinite(total):
            raise ValueError(f"{path}: total abundance is not finite")
        resolved = {}
        for identifier, value in zip(ids, values):
            transcript = model["transcript_aliases"].get(identifier, identifier)
            if transcript in resolved:
                raise ValueError(
                    f"{path}: duplicate quantification target after alias resolution: {transcript}"
                )
            resolved[transcript] = value
        missing = set(model["transcript_map"]) - set(resolved)
        if missing:
            raise ValueError(
                f"{path}: missing expected quantification targets: {sorted(missing)[:10]}"
            )
        for transcript, gene in model["transcript_map"].items():
            abundance[positions[gene], sample_index] += resolved[transcript]
        totals.append(total)
        hashes.append(dict(sample_id=row.sample_id, sha256=sha256(path)))
    metadata["assigned_abundance"] = totals
    metadata["ready"] = [(t > 0 and t >= min_reads) for t in totals]
    metadata.to_csv(output / "sample_metadata.tsv", sep="\t", index=False)
    provenance = dict(
        schema_version=1,
        batch_index=batch_index,
        analysis_index=analysis_index,
        report_sequence=report_sequence,
        model=model,
        quantifications=hashes,
        manifest_sha256=sha256(manifest),
        settings=dict(
            log_base=log_base,
            pseudocount=pseudocount,
            min_read_count=min_reads,
            padj_cutoff=cutoff,
            normalization="all Oarfish assigned abundance per million; no length correction",
        ),
        control_sample_ids=metadata.loc[control, "sample_id"].tolist(),
        software=dict(
            imodulon_analysis=__version__,
            numpy=np.__version__,
            scipy=scipy.__version__,
            pandas=pd.__version__,
        ),
    )
    write_json(output / "provenance.json", provenance)
    if not metadata.ready.all():
        write_json(
            output / "status.json",
            dict(
                status="deferred",
                reason="insufficient_assigned_abundance",
                sample_ids=metadata.loc[~metadata.ready, "sample_id"].tolist(),
                statistical_availability="unavailable",
            ),
        )
        return
    abundance = abundance / np.asarray(totals) * 1e6
    logged = np.log(abundance + pseudocount) / np.log(log_base)
    reference = logged[:, control].mean(axis=1)
    centered = logged - reference[:, None]
    activities = inverse @ centered
    raw_activities = inverse @ logged
    residual = centered - m @ activities
    if not all(
        np.isfinite(x).all() for x in (logged, activities, raw_activities, residual)
    ):
        raise ValueError(
            "Projection produced nonfinite values; check matrix scaling and transformation"
        )
    ids = metadata.sample_id.tolist()
    pd.DataFrame(
        activities,
        index=pd.Index(model["components"], name="component_id"),
        columns=ids,
    ).to_csv(output / "activities.tsv", sep="\t")
    pd.DataFrame(
        centered, index=pd.Index(model["genes"], name="model_gene_id"), columns=ids
    ).to_csv(output / "centered_expression.tsv", sep="\t")
    pd.DataFrame(
        dict(
            model_gene_id=model["genes"],
            reference_expression=reference,
            control_sample_ids=json.dumps(metadata.loc[control, "sample_id"].tolist()),
        )
    ).to_csv(output / "reference_expression.tsv", sep="\t", index=False)
    long = (
        pd.DataFrame(
            activities,
            index=pd.Index(model["components"], name="component_id"),
            columns=ids,
        )
        .reset_index()
        .melt(id_vars="component_id", var_name="sample_id", value_name="activity")
        .merge(
            metadata[["sample_id", "alias", "group", "order"]],
            on="sample_id",
            validate="many_to_one",
        )
    )
    long.to_csv(output / "activities_long.tsv", sep="\t", index=False)
    long.groupby(["component_id", "group"], sort=False).activity.agg(
        ["mean", "std", "count"]
    ).rename(columns={"std": "sd", "count": "n"}).reset_index().to_csv(
        output / "activity_summary.tsv", sep="\t", index=False
    )
    differential = differential_activity(
        raw_activities, activities, metadata, model["components"], control_group, cutoff
    )
    differential.to_csv(output / "differential_activity.tsv", sep="\t", index=False)
    error = np.sum(residual**2, axis=0)
    denominator = np.sum(centered**2, axis=0)
    ratio = np.divide(
        error, denominator, out=np.full_like(error, np.nan), where=denominator > 0
    )
    pd.DataFrame(
        dict(
            sample_id=ids,
            residual_sum_squares=error,
            centered_sum_squares=denominator,
            residual_rmse=np.sqrt(error / len(model["genes"])),
            normalized_residual=ratio,
        )
    ).to_csv(output / "projection_qc.tsv", sep="\t", index=False)
    for filename in ("gene_mapping.tsv", "component_coverage.tsv"):
        (output / filename).write_bytes((prepared / filename).read_bytes())
    tested = int(differential.status.eq("tested").sum())
    write_json(
        output / "status.json",
        dict(
            status="ready",
            tested_count=tested,
            untested_count=len(differential) - tested,
            statistical_availability="complete"
            if tested and tested == len(differential)
            else "partial"
            if tested
            else "unavailable",
            statistical_status_counts=differential.status.value_counts().to_dict(),
        ),
    )


def finite_number(value, lower, upper=None, inclusive=False):
    number = float(value)
    if (
        not math.isfinite(number)
        or (number < lower if inclusive else number <= lower)
        or (upper is not None and number > upper)
    ):
        raise argparse.ArgumentTypeError(f"Invalid numeric value: {value}")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    prep.add_argument("--matrix", required=True)
    prep.add_argument("--annotation", required=True)
    prep.add_argument("--gene-map")
    prep.add_argument(
        "--min-gene-coverage", type=lambda x: finite_number(x, 0, 1), default=1.0
    )
    prep.add_argument("--output", required=True)
    run = modes.add_parser("analyze")
    for flag in ("prepared", "manifest", "counts-dir", "output"):
        run.add_argument("--" + flag, required=True)
    run.add_argument("--log-base", type=lambda x: finite_number(x, 1), default=2.0)
    run.add_argument("--pseudocount", type=lambda x: finite_number(x, 0), default=1.0)
    run.add_argument(
        "--min-reads",
        type=lambda x: finite_number(x, 0, inclusive=True),
        default=10000.0,
    )
    run.add_argument("--cutoff", type=lambda x: finite_number(x, 0, 1), default=0.05)
    for flag in ("batch-index", "analysis-index", "report-sequence"):
        run.add_argument("--" + flag, type=int, default=0)
    args = vars(parser.parse_args())
    mode = args.pop("mode")
    if mode == "prepare":
        args["min_coverage"] = args.pop("min_gene_coverage")
        prepare(**args)
    else:
        analyze(**args)
