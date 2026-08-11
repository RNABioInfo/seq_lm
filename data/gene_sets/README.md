# Pseudomonas putida carbon-stress gene sets

These GMT collections reproduce the 12 functional gene sets in the `Gene sets`
worksheet of Data S1 from:

> Ankenbauer A, Schäfer RA, Viegas SC, et al. (2020).
> *Pseudomonas putida* KT2440 is naturally endowed to withstand industrial-scale
> stress conditions. Microbial Biotechnology 13:1145–1161.
> <https://doi.org/10.1111/1751-7915.13571>

The experiment compared cells sampled after acute glucose depletion in a plug
flow reactor with their glucose-limited stirred-tank controls. The sets should
therefore be interpreted as pathways tested in a glucose-starvation experiment,
not as a universal signature of every form of carbon limitation.

## Files

- `p_putida_kt2440_carbon_stress.gmt`: KT2440 `PP_####` locus tags.
- `p_putida_kt2440_carbon_stress_np.gmt`: KT2440 RefSeq `NP_...` accessions,
  matching the published count matrix.
- `p_putida_nbrc14164_carbon_stress.gmt`: ortholog-translated NBRC 14164
  `PP4_RS...` locus tags. In this annotation, `gene_id` and `locus_tag` are
  identical. This is the recommended, biologically interpretable collection.
- `p_putida_nbrc14164_carbon_stress_transcript_ids.gmt`: the same collection
  expressed as Oarfish `tname` values (`unassigned_transcript_...`) for tools
  that cannot resolve locus tags through a GTF.
- `p_putida_kt2440_to_nbrc14164_mapping.tsv`: auditable protein and transcript
  mapping.
- `p_putida_carbon_stress_gene_set_summary.tsv`: set sizes and mapping coverage.

The NBRC 14164 translation used reciprocal best BLASTP hits between the complete
KT2440 (`GCF_000007565.2`) and NBRC 14164 (`GCF_000412675.1`) proteomes. Hits
were required to have at least 70% amino-acid identity and at least 70% query
and subject coverage. Of 334 unique source proteins, 307 mapped unambiguously
to 305 unique NBRC 14164 loci (two pairs of KT2440 proteins converge on the
same NBRC loci). Per-set membership retention ranges from 75% to 100%; consult
the summary table before interpreting a weak or null result.

## Running `fry`

The analysis resolves GMT locus tags to the `tname` column in the Oarfish
quantification files through `--annotation`. For data quantified with the
bundled NBRC 14164 reference:

```bash
edgeR-analysis \
    --quant_manifest quant_manifest.tsv \
    --annotation genomic.oarfish.gtf \
    --gene_sets data/gene_sets/p_putida_nbrc14164_carbon_stress.gmt \
    --output_dir edgeR_results
```

For example, locus tag `PP4_RS04630` in the recommended GMT resolves through
the GTF to `transcript_id "unassigned_transcript_947"`, which is the identifier
reported in the Oarfish `tname` column. If `--annotation` cannot be supplied,
use `p_putida_nbrc14164_carbon_stress_transcript_ids.gmt` instead.

The command reports the number of GMT memberships resolved to the count
matrix and stops with an identifier-specific error if none resolve. Always
inspect `gene_set_coverage.tsv`: `count_matrix_coverage` diagnoses identifier
mapping, whereas `tested_coverage` additionally reflects expression filtering.

The analysis tests every non-control group against `control`. Override the
reference label with `--control_group`. Results include:

- `gene_set_resolution.tsv`, the member-by-member identifier mapping;
- `gene_set_coverage.tsv`, written once for the analysis;
- `edgeR_mds_data.tsv`, the sample coordinates behind the overview MDS plot;
- `edgeR_bcv_data.tsv`, the feature-level values behind the overview BCV plot;
- `group_<target>_vs_<control>/edgeR_results.tsv`;
- `group_<target>_vs_<control>/fry_results.tsv`;
- `group_<target>_vs_<control>/fry_signed_significance.png`.

The integrated `qc_report.html` adds a Gene Set Enrichment section. It shows a
signed directional-FDR summary for each contrast and an offline dropdown with
one logFC-ranked barcode and limma-style enrichment worm per tested gene set.
The static PNG is retained as a standalone edgeR output.

Use the directional `FDR` column for coordinated up- or downregulation. The
`FDR.Mixed` column tests whether genes in a set change in either direction and
can be significant even when the directional test is not. Enrichment indicates
association with the expression contrast; it does not establish pathway
activity or causation.

The source article and supplementary data are distributed under
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
