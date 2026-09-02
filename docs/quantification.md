# Quantification and transcript-biotype QC

Oarfish quantification is enabled automatically when both
`--reference_genome` and `--reference_annotation` are supplied. The two
reference inputs form a pair: supplying only one is an error. Quantification
does not depend on `--differential_expression`; with differential expression
disabled, the workflow still publishes per-sample quantifications, persists
and restores final-sample checkpoints, and adds a Transcript biotypes tab to
the Quality Control report.

For each complete live batch, the workflow name-collates each new BAM chunk,
assembles each sample's cumulative alignment stream, and runs Oarfish once for
that cumulative sample snapshot. The newest quantification for every sample is
then emitted as an ordered, complete experiment batch. Differential expression,
when enabled, consumes these completed batches rather than launching a second
Oarfish run.

Quantifications are published under:

```text
<out_dir>/<group>/<alias>/quantification/
```

## Transcript-biotype classification

The annotation is parsed once per workflow run. Oarfish `tname` identifiers are
matched first to exact `transcript_id` values. GFF3 `ID` aliases and an
unambiguous identifier with a standard namespace prefix removed (for example,
`transcript:tx1` to `tx1`) are also supported. GFF3 attribute values are URL
decoded.

Biotype resolution uses the first available annotation level in this order:

1. transcript `transcript_biotype`, `transcript_type`, then `biotype`;
2. transcript-linked gene `gene_biotype`, `gene_type`, then `biotype`.

Missing, unmatched, or conflicting classifications are `Unknown`; the workflow
does not infer a biotype from expression, feature names, or sequence content.
Nonempty labels are normalized to this fixed report order:

1. `Protein-coding`
2. `rRNA`
3. `tRNA`
4. `lncRNA`
5. `Other ncRNA`
6. `Pseudogene`
7. `Other`
8. `Unknown`

`Other ncRNA` includes recognized classes such as miRNA, snRNA, snoRNA, SRP
RNA, RNase P RNA, and tmRNA. Unrecognized nonempty labels are `Other`.

For every sample, the workflow sums Oarfish's EM-estimated `num_reads` by the
eight canonical classes and divides by the sum across all classes, including
`Unknown`. These are fractions among transcript-assigned Oarfish abundance—not
fractions of all sequenced or mapped reads. A zero-total sample has zero for
every fraction and is explicitly marked as having no assigned reads. Malformed,
negative, or non-finite Oarfish abundance values stop the analysis.

The report input TSV has the stable contract:

```text
name	group	biotype	num_reads	fraction
```

All eight rows are emitted for every sample, including zero-valued groups. The
optional `workflow-glue qc_report --transcript-biotypes` argument consumes this
table; omitting it preserves the existing report layout.

## Prokaryotic annotations

NCBI prokaryotic GTF annotations commonly lack the transcript and exon records
required by Oarfish for protein-coding targets. Convert those annotations with
`bin/oarfish-gtf-convert` as described in the
[differential-expression documentation](differential_expression.md#ncbi-prokaryotic-annotations).
The converter retains declared RNA models and carries declared gene/transcript
biotypes into the Oarfish-compatible representation. It does not infer operons
or transcript boundaries.
