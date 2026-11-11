# Quantification and transcript biotypes

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

Oarfish estimates transcript abundance from aligned reads. Enable
quantification by setting both reference paths in the `params` block of
the run configuration.

```groovy
reference_genome = "/data/reference.fa"
reference_annotation = "/data/transcripts.gtf"
```

Use the reference FASTA from the alignment step and a matching GTF or GFF3
annotation with transcript and exon records. To run only quality control and
quantification, set `differential_expression`, `gene_set_enrichment`, and
`ica_analysis` to `false`.

## Strand specificity

Set `strand_specificity` according to the library preparation. The selected
setting applies to all samples in the run.

| `strand_specificity` | Read orientation relative to the transcript |
| --- | --- |
| `forward` (default) | Same orientation |
| `reverse` | Opposite orientation |
| `unstranded` | Both orientations |

## Prokaryotic annotations

Convert NCBI GTF annotations that provide coding sequence records without
transcript and exon records before running seq_lm.

```bash
bin/oarfish-gtf-convert genomic.gtf genomic.oarfish.gtf
```

Set `reference_annotation` to `genomic.oarfish.gtf`. The converter preserves
existing RNA models and creates one transcript spanning each protein-coding
gene. Use a curated annotation for analyses that require operon boundaries
or specific transcript isoforms.

## Results

During live analysis, each abundance estimate uses all reads processed for
that sample to date. The final estimates are available in
`<out_dir>/<group>/<alias>/quantification/final.quant`.

In the report, **Quality Control → Transcript Biotypes** shows the proportion
of Oarfish-assigned reads in each transcript biotype. The reference annotation
provides the labels. Reads assigned to unclassified transcripts count toward
the total and appear under **Unknown**.
