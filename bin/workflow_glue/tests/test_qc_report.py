"""Test qc_report.py."""

import gzip
import json

import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from workflow_glue import qc_report  # noqa: E402
from workflow_glue.qc_report_types import histogram_plots, kde_plots  # noqa: E402
from workflow_glue.transcript_biotypes import BIOTYPE_ORDER  # noqa: E402


def test_write_report_uses_live_shell_and_versioned_snapshot(tmp_path):
    """Live reports refresh their frame and preserve active tab selections."""

    class Report:
        def write(self, path):
            path = tmp_path / path if not hasattr(path, "write_text") else path
            path.write_text("<html><body>snapshot</body></html>")

    report_path = tmp_path / "qc_report.html"
    qc_report.write_report(Report(), report_path, "batch 2", 5)

    shell = report_path.read_text()
    state = json.loads((tmp_path / "qc_report_state.json").read_text())
    snapshot = tmp_path / state["snapshot"]
    assert snapshot.name == "qc_report_snapshot_batch_2.html"
    snapshot_html = snapshot.read_text()
    assert "snapshot" in snapshot_html
    assert 'id="seq-lm-report-navigation"' in snapshot_html
    assert state["latest_batch"] == "batch 2"
    assert state["snapshot_bytes"] == snapshot.stat().st_size
    assert '<meta http-equiv="refresh"' not in shell
    assert "captureViewState" in shell
    assert "restoreViewState" in shell
    assert "seqLmTabKey" in shell
    assert "Object.fromEntries" in shell
    assert "frameReady" in shell
    assert "stateIsOlder" in shell
    assert "fetch(requestUrl" in shell
    assert 'method: "HEAD"' in shell
    assert "Waiting for update." in shell
    assert "Waiting for the next complete report update" not in shell


def test_write_report_rebrands_labs_header(tmp_path):
    """Generated reports use the project logo, color, and subheadline."""

    class Report:
        def write(self, path):
            path.write_text(
                "<html><head></head><body><header>"
                '<nav class="fixed-top bg-dark">'
                '<a href="https://labs.epi2me.io/">'
                '<div alt="EPI2ME Labs Logo"><svg></svg></div>'
                "</a></nav>"
                '<p class="py-2 fs-5">Results generated through the qc_report '
                "workflow provided by Oxford Nanopore Technologies.</p>"
                "</header></body></html>"
            )

    report_path = tmp_path / "qc_report.html"
    qc_report.write_report(Report(), report_path, "2", 0)

    html = report_path.read_text()
    assert "Waiting for update." not in html
    assert "https://labs.epi2me.io/" not in html
    assert "EPI2ME Labs Logo" not in html
    assert 'alt="RNA BioInfo AUCG logo"' in html
    assert "data:image/png;base64," in html
    assert "#004191" in html
    assert (
        "Long-read RNA sequencing quality control, differential expression, "
        "and gene set enrichment."
    ) in html
    assert "Oxford Nanopore Technologies" not in html


@pytest.mark.parametrize(
    "notice",
    [
        "For DEA, the required read depth is not yet satisfied.",
        "For DEA, the sample quantifications contain no matching feature IDs.",
    ],
)
def test_write_report_adds_dea_readiness_notice_to_subtitle(tmp_path, notice):
    """A failed DEA precondition is visible without adding a DEA tab."""

    class Report:
        def write(self, path):
            path.write_text(
                '<html><head></head><body><header><h1>Report</h1>'
                '<p class="py-2 fs-5">Workflow results.</p>'
                '</header><main>Quality Control</main></body></html>'
            )

    report_path = tmp_path / "qc_report.html"
    qc_report.write_report(
        Report(),
        report_path,
        "2",
        0,
        notice,
    )

    html = report_path.read_text()
    assert "Workflow results." not in html
    assert (
        "Long-read RNA sequencing quality control, differential expression, "
        "and gene set enrichment."
    ) in html
    assert notice in html
    assert 'class="seq-lm-dea-readiness-notice"' in html


def write_flagstat(path):
    """Write enough samtools flagstat rows for the QC report parser."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "100\t0\ttotal (QC-passed reads + QC-failed reads)\n"
        "80\t0\tprimary\n"
        "15\t0\tsecondary\n"
        "5\t0\tsupplementary\n"
        "0\t0\tduplicates\n"
        "0\t0\tprimary duplicates\n"
        "60\t0\tmapped\n"
        "60.00%\tN/A\tmapped %\n"
        "55\t0\tprimary mapped\n"
    )


def write_nanoplot(path, rows):
    """Write a gzipped NanoPlot-data style TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(
            "readIDs\tquals\taligned_quals\tlengths\taligned_lengths\tmapQ\t"
            "percentIdentity\n"
        )
        for row in rows:
            handle.write(
                "{readIDs}\t{quals}\t{aligned_quals}\t{lengths}\t"
                "{aligned_lengths}\t{mapQ}\t{percentIdentity}\n".format(**row)
            )


def nanoplot_df(rows):
    """Build a NanoPlot-like DataFrame for focused report tests."""
    return pd.DataFrame(rows)


def write_transcript_biotypes(path):
    """Write a complete fixed-category composition table for two samples."""
    rows = []
    for name, group, protein_fraction in (
        ("control_1", "control", 0.75),
        ("treatment_1", "time_point_1", 0.5),
    ):
        for biotype in BIOTYPE_ORDER:
            fraction = protein_fraction if biotype == "Protein-coding" else 0.0
            if biotype == "Unknown":
                fraction = 1.0 - protein_fraction
            rows.append(
                {
                    "name": name,
                    "group": group,
                    "biotype": biotype,
                    "num_reads": fraction * 100,
                    "fraction": fraction,
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_differential_results(path):
    """Write a small report-ready edgeR result tree."""
    path.mkdir(parents=True)
    (path / "feature_counts.tsv").write_text(
        "feature_id\tcontrol_1\ttreatment_1\n"
        "gene_1\t2\t32\n"
        "gene_2\t40\t4\n"
        "gene_3\t10\t12\n"
        "gene_4\t8\t9\n"
    )
    (path / "edgeR_bcv_data.tsv").write_text(
        "feature_id\taverage_log_cpm\ttagwise_dispersion\ttagwise_bcv\t"
        "trended_dispersion\ttrended_bcv\tcommon_dispersion\tcommon_bcv\n"
        "gene_1\t3\t0.16\t0.4\t0.1225\t0.35\t0.09\t0.3\n"
        "gene_2\t4\t0.09\t0.3\t0.09\t0.3\t0.09\t0.3\n"
        "gene_3\t5\t0.04\t0.2\t0.0625\t0.25\t0.09\t0.3\n"
        "gene_4\t6\t0.01\t0.1\t0.04\t0.2\t0.09\t0.3\n"
    )
    (path / "edgeR_mds_data.tsv").write_text(
        "sample\tgroup\tdimension_1\tdimension_2\tdimension_1_variance\t"
        "dimension_2_variance\taxis_label\ttop_features\tgene_selection\n"
        "control_1\tcontrol\t-1\t0.5\t0.7\t0.2\tLeading logFC dim\t4\tpairwise\n"
        "treatment_1\ttime_point_1\t1\t-0.5\t0.7\t0.2\tLeading logFC dim\t4\tpairwise\n"
    )
    contrast = path / "group_time_point_1_vs_control"
    contrast.mkdir()
    (contrast / "edgeR_results.tsv").write_text(
        "feature_id\tlogFC\tlogCPM\tPValue\tFDR\tgene\n"
        "gene_1\t2.5\t5.0\t0.0001\t0.001\tup_gene\n"
        "gene_2\t-2.2\t5.5\t0.0002\t0.002\tdown_gene\n"
        "gene_3\t0.3\t4.0\t0.5\t0.8\tneutral_gene\n"
        "gene_4\t0.0\t3.0\t1.0\t1.0\tconstant_gene\n"
    )
    (contrast / "fry_results.tsv").write_text(
        "gene_set\tNGenes\tDirection\tPValue\tFDR\tPValue.Mixed\t"
        "FDR.Mixed\tdescription\tgmt_members\tmatched_gmt_members\t"
        "count_matrix_members\ttested_members\ttested_gmt_members\t"
        "count_matrix_coverage\ttested_coverage\n"
        "carbon_up\t2\tUp\t0.001\t0.01\t0.002\t0.02\t"
        "Carbon response\t2\t2\t2\t2\t2\t1.0\t1.0\n"
        "carbon_mixed\t2\tDown\t0.8\t0.8\t0.001\t0.01\t"
        "Mixed response\t2\t2\t2\t2\t2\t1.0\t1.0\n"
    )
    (path / "gsva_scores_long.tsv").write_text(
        "gene_set\tdescription\tn_genes\tsample\tgroup\tscore\n"
        "carbon_up\tCarbon response\t2\tcontrol_1\tcontrol\t-0.6\n"
        "carbon_up\tCarbon response\t2\ttreatment_1\ttime_point_1\t0.7\n"
        "carbon_mixed\tMixed response\t2\tcontrol_1\tcontrol\t0.4\n"
        "carbon_mixed\tMixed response\t2\ttreatment_1\ttime_point_1\t-0.3\n"
    )
    (path / "gsva_gene_set_coverage.tsv").write_text(
        "gene_set\tdescription\tresolved_members\tretained_members\t"
        "variable_members\tscored_members\tstatus\n"
        "carbon_up\tCarbon response\t2\t2\t2\t2\tscored\n"
        "carbon_mixed\tMixed response\t2\t2\t2\t2\tscored\n"
    )
    (contrast / "gsva_limma_results.tsv").write_text(
        "gene_set\tdescription\tn_genes\ttarget_group\tcontrol_group\t"
        "effect_size\taverage_score\tt_statistic\tp_value\t"
        "adjusted_p_value\tlog_odds\n"
        "carbon_up\tCarbon response\t2\ttime_point_1\tcontrol\t"
        "1.3\t0.05\t5.0\t0.001\t0.002\t2.0\n"
        "carbon_mixed\tMixed response\t2\ttime_point_1\tcontrol\t"
        "-0.7\t0.05\t-3.0\t0.01\t0.01\t1.0\n"
    )
    (path / "gene_set_resolution.tsv").write_text(
        "gene_set\tfeature_id\n"
        "carbon_up\tgene_1\n"
        "carbon_up\tgene_4\n"
        "carbon_mixed\tgene_2\n"
        "carbon_mixed\tgene_3\n"
    )


def write_stability_results(path):
    """Write the matching sample-level stability audit."""
    path.write_text(
        "analysis_index\tbatch_index\tgroup\tsample\tbam_dir\t"
        "effectively_live\trequired_contrasts\tconsecutive_stable_batches\t"
        "eligible\tnewly_eligible\tbehavior\taction_result\n"
        "3\t2\tcontrol\tcontrol_1\t/data/control_1\ttrue\t"
        "group_time_point_1_vs_control\t2\tfalse\tfalse\tlog\tnone\n"
        "3\t2\ttime_point_1\ttreatment_1\t/data/treatment_1\ttrue\t"
        "group_time_point_1_vs_control\t3\ttrue\ttrue\tlog\tlogged\n"
    )


def test_create_hist_plot_uses_length_values_only(monkeypatch):
    """Length histogram passes one cleaned numeric column to ezCharts."""
    calls = {}

    class Axis:
        axis_label = None

    class Figure:
        xaxis = Axis()
        yaxis = Axis()

    class Plot:
        title = None
        _fig = Figure()

    def fake_histplot(data, bins):
        calls["data"] = data
        calls["bins"] = bins
        return Plot()

    monkeypatch.setattr(histogram_plots.ezc, "histplot", fake_histplot)

    plot = histogram_plots.create_hist_plot(
        nanoplot_df(
            [
                {"lengths": "100", "mapQ": "30"},
                {"lengths": "not-a-number", "mapQ": "40"},
                {"lengths": "250", "mapQ": "50"},
            ]
        ),
        x_column="lengths",
        title="Read length histogram",
        x_axis_label="Read length",
        y_axis_label="Number of reads",
        bins=25,
    )

    assert calls["data"].tolist() == [100, 250]
    assert calls["bins"] == 25
    assert plot.title == {"text": "Read length histogram"}
    assert plot._fig.xaxis.axis_label == "Read length"
    assert plot._fig.yaxis.axis_label == "Number of reads"


def test_qc_report_writes_html(tmp_path):
    """Report lists the current QC result samples and QC input directories."""
    write_flagstat(
        tmp_path / "qc_results/flagstat/control/control_1/flagstat_data_chunk_0.tsv"
    )
    write_flagstat(
        tmp_path
        / "qc_results/flagstat/time_point_1/treatment_1/flagstat_data_chunk_0.tsv"
    )
    write_nanoplot(
        tmp_path / "qc_results/nanoplot/control/control_1/nanoplot_data_chunk_0.tsv.gz",
        [
            {
                "readIDs": "read_1",
                "quals": 10,
                "aligned_quals": 9,
                "lengths": 100,
                "aligned_lengths": 90,
                "mapQ": 30,
                "percentIdentity": 98.0,
            },
            {
                "readIDs": "read_2",
                "quals": 20,
                "aligned_quals": 18,
                "lengths": 200,
                "aligned_lengths": 180,
                "mapQ": 40,
                "percentIdentity": 99.0,
            },
        ],
    )
    write_nanoplot(
        tmp_path / "qc_results/nanoplot/control/control_1/nanoplot_data_chunk_1.tsv.gz",
        [
            {
                "readIDs": "read_3",
                "quals": 30,
                "aligned_quals": 27,
                "lengths": 300,
                "aligned_lengths": 270,
                "mapQ": 50,
                "percentIdentity": 97.0,
            },
        ],
    )
    write_nanoplot(
        tmp_path
        / "qc_results/nanoplot/time_point_1/treatment_1/nanoplot_data_chunk_0.tsv.gz",
        [
            {
                "readIDs": "read_4",
                "quals": 15,
                "aligned_quals": 12,
                "lengths": 150,
                "aligned_lengths": 140,
                "mapQ": 35,
                "percentIdentity": 96.0,
            },
            {
                "readIDs": "read_5",
                "quals": 25,
                "aligned_quals": 22,
                "lengths": 250,
                "aligned_lengths": 240,
                "mapQ": 45,
                "percentIdentity": 95.0,
            },
        ],
    )
    samples = tmp_path / "qc_report_samples.tsv"
    samples.write_text(
        "name\tgroup\torder\tchunks_seen\tlatest_batch_index\tqc_dir\n"
        "control_1\tcontrol\t0\t2\t2\tqc_results\n"
        "treatment_1\ttime_point_1\t15\t1\t1\tqc_results\n"
    )
    differential_results = tmp_path / "differential_results"
    write_differential_results(differential_results)
    stability_results = tmp_path / "sample_stability.tsv"
    write_stability_results(stability_results)
    transcript_biotypes = tmp_path / "transcript_biotypes.tsv"
    write_transcript_biotypes(transcript_biotypes)

    versions = tmp_path / "versions"
    versions.mkdir()
    (versions / "versions.txt").write_text("qc_report,workflow\n")

    params = tmp_path / "params.json"
    params.write_text(json.dumps({"qc_report": True}))

    report = tmp_path / "qc_report.html"
    args = qc_report.argparser().parse_args(
        [
            str(report),
            "--samples",
            str(samples),
            "--versions",
            str(versions),
            "--params",
            str(params),
            "--latest-batch",
            "2",
            "--differential-results",
            str(differential_results),
            "--stability-results",
            str(stability_results),
            "--stability-behavior",
            "log",
            "--transcript-biotypes",
            str(transcript_biotypes),
            "--gene-set-enrichment",
            "--temporal-analysis",
            "--refresh-seconds",
            "0",
        ]
    )

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        qc_report.main(args)

    html = report.read_text()
    assert "seq_lm analysis report" in html
    assert (
        html.find("Quality Control")
        < html.find("Differential Analysis")
        < html.find("Gene Set Enrichment")
        < html.find("Temporal Analysis")
    )
    assert "control/control_1" in html
    assert "time_point_1/treatment_1" in html
    assert "qc_results" in html
    assert "Number of mapped reads" in html
    assert "Number of bases" in html
    assert "Number of aligned bases" in html
    assert "Median read length" in html
    assert "Mean MapQ" in html
    assert "600" in html
    assert "540" in html
    assert "40.00" in html
    assert "Plotly.newPlot" not in html
    assert "window.PlotlyConfig" not in html
    assert "updatemenus" not in html
    assert "'type': 'sankey'" in html
    assert "Transcript biotype composition" in html
    assert "Unknown denotes targets without one unambiguous annotation biotype." in html
    assert "Fractions are EM-estimated Oarfish abundance" not in html
    assert "Protein-coding" in html
    assert "Unknown" in html
    assert (
        html.find(">Read flow<")
        < html.find(">Transcript biotypes<")
        < html.find(">Metrics<")
    )
    assert "Read length vs Read quality" in html
    assert "PCA of log2(CPM + 1)" in html
    assert "edgeR MDS of leading logFC" in html
    assert "edgeR biological coefficient of variation (BCV)" in html
    assert "logFC vs logCPM" in html
    assert "Volcano plot" in html
    assert "Top differential genes" in html
    assert "Result Stability" in html
    assert "#Stable consec. batches" in html
    assert "fry signed directional significance" in html
    assert "Gene-set barcode and enrichment worm" in html
    assert "GSVA scores across samples" in html
    assert "Raw GSVA scores" in html
    assert "GSVA limma volcano" in html
    assert "GSVA score difference" in html
    assert "GSVA limma across contrasts" in html
    assert "Differential GSVA scores" in html
    assert "GSVA score over time" in html
    assert "Gene expression over time" in html
    assert "Gene z-score" in html
    assert "Descriptive only. Heatmap colors are gene-wise z-scores." in html
    assert "Time (min)" in html
    assert "Carbon response" in html
    assert "Mixed response" in html
    assert "Relative local enrichment" in html
    assert "not a formal GSEA running-sum score" not in html
    assert "fry enrichment is an expression association" not in html
    assert "Gene-set details" in html
    assert "BEGIN bokeh.min.js" in html
    assert "BEGIN bokeh-widgets.min.js" in html
    assert html.index('t.version="5.3.3"') < html.index(
        'id="seq-lm-report-navigation"'
    )
    assert html.index("BEGIN bokeh.min.js") < html.index("BEGIN bokeh-widgets.min.js")
    assert "time_point_1 vs control" in html
    assert "control/control_1" in html
    assert "time_point_1/treatment_1" in html
    assert "dropdown-menu" in html
    assert "seq-lm-primary-tablist" in html
    assert (
        "#pills-tab.seq-lm-primary-tablist > .nav-item > .nav-link.active"
        in html
    )
    assert "shown.bs.tab" in html
    assert "resize_layout" in html
    assert "getInstanceByDom" in html

    differential_only_report = tmp_path / "qc_report_differential_only.html"
    differential_only_args = qc_report.argparser().parse_args(
        [
            str(differential_only_report),
            "--samples",
            str(samples),
            "--versions",
            str(versions),
            "--params",
            str(params),
            "--differential-results",
            str(differential_results),
            "--refresh-seconds",
            "0",
        ]
    )
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        qc_report.main(differential_only_args)
    differential_only_html = differential_only_report.read_text()
    assert ">Differential Analysis<" in differential_only_html
    assert ">Gene Set Enrichment<" not in differential_only_html
    assert ">Temporal Analysis<" not in differential_only_html

    qc_only_report = tmp_path / "qc_report_qc_only.html"
    qc_only_args = qc_report.argparser().parse_args(
        [
            str(qc_only_report),
            "--samples",
            str(samples),
            "--versions",
            str(versions),
            "--params",
            str(params),
            "--dea-readiness-notice",
            "For DEA, the required read depth is not yet satisfied.",
            "--refresh-seconds",
            "0",
        ]
    )
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        qc_report.main(qc_only_args)
    qc_only_html = qc_only_report.read_text()
    assert ">Quality Control<" in qc_only_html
    assert ">Differential Analysis<" not in qc_only_html
    assert ">Gene Set Enrichment<" not in qc_only_html
    assert ">Temporal Analysis<" not in qc_only_html
    assert "Transcript biotype composition" not in qc_only_html
    assert "For DEA, the required read depth is not yet satisfied." in qc_only_html

    biotype_only_report = tmp_path / "qc_report_biotype_only.html"
    biotype_only_args = qc_report.argparser().parse_args(
        [
            str(biotype_only_report),
            "--samples",
            str(samples),
            "--versions",
            str(versions),
            "--params",
            str(params),
            "--transcript-biotypes",
            str(transcript_biotypes),
            "--refresh-seconds",
            "0",
        ]
    )
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        qc_report.main(biotype_only_args)
    biotype_only_html = biotype_only_report.read_text()
    assert ">Quality Control<" in biotype_only_html
    assert ">Differential Analysis<" not in biotype_only_html
    assert "Transcript biotype composition" in biotype_only_html
    assert ">Transcript biotypes<" in biotype_only_html


def test_gene_set_report_requires_differential_results():
    """Gene-set tabs cannot be requested without differential results."""
    args = qc_report.argparser().parse_args(
        [
            "report.html",
            "--samples",
            "samples.tsv",
            "--versions",
            "versions",
            "--params",
            "params.json",
            "--gene-set-enrichment",
        ]
    )

    with pytest.raises(
        ValueError,
        match="--gene-set-enrichment requires --differential-results",
    ):
        qc_report.main(args)


def test_temporal_report_requires_gene_set_enrichment():
    """Temporal plots cannot be requested without their GSVA inputs."""
    args = qc_report.argparser().parse_args(
        [
            "report.html",
            "--samples",
            "samples.tsv",
            "--versions",
            "versions",
            "--params",
            "params.json",
            "--temporal-analysis",
        ]
    )

    with pytest.raises(
        ValueError,
        match="--temporal-analysis requires --gene-set-enrichment",
    ):
        qc_report.main(args)


def test_create_read_fate_sankey_html_returns_embeddable_fragment():
    """Read fate Sankey is returned as native ECharts options JSON."""
    flagstat = qc_report.FlagstatResult(
        total_reads=100,
        primary_reads=80,
        secondary_reads=15,
        supplementary_reads=5,
        total_mapped=72,
        primary_mapped=60,
    )

    html = qc_report.create_read_fate_sankey_html(flagstat)
    options = json.loads(html)

    assert options["series"][0]["type"] == "sankey"
    assert "Plotly.newPlot" not in html
    assert "Total entries" in html
    assert "Primary entries" in html
    assert "Secondary entries" in html
    assert "Primary mapped" in html
    assert "Primary unmapped" in html
    assert '"x":' not in html
    assert '"y":' not in html
    assert "Supplementary" not in html
    assert "Non-primary" not in html
    assert "<html" not in html


def test_create_sample_read_fate_sankey_html_adds_sample_dropdown():
    """Multi-sample Sankey options include all sample labels."""
    sample_results = [
        qc_report.SampleQCResult(
            "rep_1",
            "control",
            qc_report.FlagstatResult(100, 80, 15, 5, 60, 55),
            nanoplot_df([]),
        ),
        qc_report.SampleQCResult(
            "rep_1",
            "time_point_1",
            qc_report.FlagstatResult(120, 90, 20, 10, 70, 65),
            nanoplot_df([]),
        ),
    ]

    html = qc_report.create_sample_read_fate_sankey_html(sample_results)
    payload = json.loads(html)

    assert payload[0]["plot"]["series"][0]["type"] == "sankey"
    assert payload[1]["plot"]["series"][0]["type"] == "sankey"
    assert "Plotly.newPlot" not in html
    assert "updatemenus" not in html
    assert "control/rep_1" in html
    assert "time_point_1/rep_1" in html
    assert "Primary entries" in html
    assert "Secondary entries" in html
    assert '"x":' not in html
    assert '"y":' not in html
    assert "Supplementary" not in html


def test_create_read_length_quality_kde_html_handles_single_read():
    """Sparse live data renders an empty 2D KDE plot instead of raising."""
    nanoplot = nanoplot_df(
        [
            {
                "readIDs": "read_1",
                "quals": "17",
                "aligned_quals": "16",
                "lengths": "120",
                "aligned_lengths": "115",
                "mapQ": "30",
                "percentIdentity": "98.0",
            }
        ]
    )
    html = qc_report.create_read_length_quality_kde_html(nanoplot)

    assert "Read length vs quality 2D KDE (insufficient data)" in html
    assert "Plotly.newPlot" not in html
    assert "scatter" not in html


def test_kde_values_are_deterministically_limited_to_20_000_rows():
    """Large KDE inputs use the same bounded sample on every invocation."""
    values = kde_plots.np.column_stack(
        (
            kde_plots.np.arange(25_000, dtype=float),
            kde_plots.np.arange(25_000, dtype=float) % 97,
        )
    )

    first = kde_plots._sample_kde_values(values)
    second = kde_plots._sample_kde_values(values)

    assert len(first) == 20_000
    assert kde_plots.np.array_equal(first, second)
    assert len(kde_plots._sample_kde_values(values[:100])) == 100


def test_kde_default_grid_size_is_100():
    """The default KDE evaluation grid is bounded to 100 by 100 cells."""
    assert (
        kde_plots.create_2d_kde_plot.__defaults__[-1] == kde_plots.KDE_GRID_SIZE == 100
    )


def test_create_read_length_quality_kde_html_uses_2d_contour_plot():
    """2D KDE uses a Bokeh contour plot for valid numeric pairs."""
    html = qc_report.create_read_length_quality_kde_html(
        nanoplot_df(
            [
                {"lengths": "100", "quals": "10"},
                {"lengths": "not-a-number", "quals": "20"},
                {"lengths": "300", "quals": "30"},
                {"lengths": "450", "quals": "40"},
                {"lengths": "700", "quals": "26"},
            ]
        )
    )

    assert "Read length vs quality 2D KDE" in html
    assert "insufficient data" not in html
    assert "ContourRenderer" in html
    assert "Plotly.newPlot" not in html
    assert "scatter" not in html


def test_create_2d_kde_html_accepts_column_names():
    """Generic 2D KDE helper accepts caller-selected x/y columns."""
    html = qc_report.create_2d_kde_html(
        nanoplot_df(
            [
                {"aligned_lengths": "100", "mapQ": "10"},
                {"aligned_lengths": "300", "mapQ": "30"},
                {"aligned_lengths": "450", "mapQ": "40"},
                {"aligned_lengths": "700", "mapQ": "26"},
            ]
        ),
        x_column="aligned_lengths",
        y_column="mapQ",
        title="Aligned read length vs MapQ 2D KDE",
        x_axis_label="Aligned read length",
        y_axis_label="MapQ",
    )

    assert "Aligned read length vs MapQ 2D KDE" in html
    assert "Aligned read length" in html
    assert "MapQ" in html
    assert "ContourRenderer" in html
    assert "insufficient data" not in html


def test_create_nanoplot_metrics_table_summarizes_samples():
    """NanoPlot metrics are summarized with one column per sample."""
    sample_results = [
        qc_report.SampleQCResult(
            "rep_1",
            "control",
            qc_report.FlagstatResult(100, 80, 15, 5, 60, 55),
            nanoplot_df(
                [
                    {
                        "readIDs": "read_1",
                        "quals": 10,
                        "aligned_quals": 9,
                        "lengths": 100,
                        "aligned_lengths": 90,
                        "mapQ": 30,
                        "percentIdentity": 98.0,
                    },
                    {
                        "readIDs": "read_2",
                        "quals": 20,
                        "aligned_quals": 18,
                        "lengths": 200,
                        "aligned_lengths": 180,
                        "mapQ": 40,
                        "percentIdentity": 99.0,
                    },
                    {
                        "readIDs": "read_3",
                        "quals": 30,
                        "aligned_quals": 27,
                        "lengths": 300,
                        "aligned_lengths": 270,
                        "mapQ": 50,
                        "percentIdentity": 97.0,
                    },
                ]
            ),
        ),
        qc_report.SampleQCResult(
            "rep_1",
            "time_point_1",
            qc_report.FlagstatResult(120, 90, 20, 10, 70, 65),
            nanoplot_df(
                [
                    {
                        "readIDs": "read_4",
                        "quals": 15,
                        "aligned_quals": 12,
                        "lengths": 150,
                        "aligned_lengths": 140,
                        "mapQ": 35,
                        "percentIdentity": 96.0,
                    },
                    {
                        "readIDs": "read_5",
                        "quals": 25,
                        "aligned_quals": 22,
                        "lengths": 250,
                        "aligned_lengths": 240,
                        "mapQ": 45,
                        "percentIdentity": 95.0,
                    },
                ]
            ),
        ),
    ]

    metrics = qc_report.create_nanoplot_metrics_table(sample_results)
    values = metrics.set_index("Metric")

    assert list(metrics.columns) == [
        "Metric",
        "control/rep_1",
        "time_point_1/rep_1",
    ]
    assert values.loc["Number of mapped reads", "control/rep_1"] == "3"
    assert values.loc["Number of bases", "control/rep_1"] == "600"
    assert values.loc["Number of aligned bases", "control/rep_1"] == "540"
    assert values.loc["Median read length", "control/rep_1"] == "200.00"
    assert values.loc["Mean read length", "control/rep_1"] == "200.00"
    assert values.loc["Median read quals", "control/rep_1"] == "20.00"
    assert values.loc["Mean read quals", "control/rep_1"] == "20.00"
    assert values.loc["Median MapQ", "control/rep_1"] == "40.00"
    assert values.loc["Mean MapQ", "control/rep_1"] == "40.00"
    assert values.loc["Number of mapped reads", "time_point_1/rep_1"] == "2"
    assert values.loc["Number of bases", "time_point_1/rep_1"] == "400"
    assert values.loc["Number of aligned bases", "time_point_1/rep_1"] == "380"
