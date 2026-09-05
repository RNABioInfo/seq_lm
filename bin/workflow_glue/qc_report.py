"""Create the live QC, differential-analysis, and enrichment report."""

from html import escape
import json
from pathlib import Path
import re

from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable

from .report_compat import apply_report_branding, labs_report
from .util import get_named_logger, wf_parser  # noqa: ABS101

from .qc_report_types.differential_plots import (
    add_differential_analysis,
    load_differential_results,
    load_stability_results,
)
from .qc_report_types.gene_set_plots import (
    add_gene_set_enrichment,
    load_gene_set_results,
)
from .qc_report_types.gsva_plots import (
    add_gsva_differential,
    add_gsva_scores,
    load_gsva_results,
)
from .qc_report_types.temporal_plots import (
    add_temporal_analysis,
    load_temporal_results,
)
from .qc_report_types.imodulon_plots import (
    add_imodulon_analysis,
    load_imodulon_results,
)
from .qc_report_types.sankey_plot import (
    add_sample_read_fate_sankeys,
    create_read_fate_sankey_html,
    create_sample_read_fate_sankey_html,
)
from .qc_report_types.parse_inputs import load_qc_samples
from .qc_report_types.kde_plots import add_sample_2d_kdes, create_2d_kde_html
from .qc_report_types.base_metrics import create_nanoplot_metrics_table
from .qc_report_types.histogram_plots import add_sample_hists
from .qc_report_types.biotype_plot import add_transcript_biotype_composition
from .qc_report_types.result_types import FlagstatResult, SampleQCResult

__all__ = [
    "FlagstatResult",
    "SampleQCResult",
    "create_2d_kde_html",
    "create_nanoplot_metrics_table",
    "create_read_fate_sankey_html",
    "create_read_length_quality_kde_html",
    "create_sample_read_fate_sankey_html",
]


def create_read_length_quality_kde_html(nanoplot):
    """Compatibility wrapper for the focused read-length/quality probe."""
    return create_2d_kde_html(
        nanoplot,
        x_column="lengths",
        y_column="quals",
        title="Read length vs quality 2D KDE",
        x_axis_label="Read length",
        y_axis_label="Read quality",
    )


def _next_report_revision(report_path):
    """Allocate beyond published or interrupted local report snapshots."""
    report_path = Path(report_path)
    pattern = re.escape(report_path.stem) + r"_snapshot_revision_(\d+)\.html"
    revisions = [
        int(match.group(1))
        for path in report_path.parent.glob(f"{report_path.stem}_snapshot_revision_*.html")
        if (match := re.fullmatch(pattern, path.name))
    ]
    state_path = report_path.with_name(f"{report_path.stem}_state.json")
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if "report_revision" in state:
            revisions.append(int(state["report_revision"]))
    return max(revisions, default=-1) + 1


def _live_report_paths(report_path, report_revision):
    """Return stable-shell, state, and immutable snapshot paths."""
    report_path = Path(report_path)
    snapshot_path = report_path.with_name(
        f"{report_path.stem}_snapshot_revision_{report_revision}{report_path.suffix}"
    )
    state_path = report_path.with_name(f"{report_path.stem}_state.json")
    return report_path, state_path, snapshot_path


def _live_report_shell(
    state_name,
    snapshot_name,
    report_revision,
    refresh_seconds,
):
    """Return a stable page shell that refreshes only its report frame."""
    refresh_milliseconds = max(int(refresh_seconds), 1) * 1000
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape("seq_lm analysis report")}</title>
  <style>
    html, body, #seq-lm-report-frame {{
      border: 0;
      height: 100%;
      margin: 0;
      padding: 0;
      width: 100%;
    }}
    #seq-lm-live-status {{
      background: #fff3cd;
      border: 1px solid #ffecb5;
      border-radius: 0.375rem;
      bottom: 1rem;
      color: #664d03;
      display: none;
      font: 14px/1.4 system-ui, sans-serif;
      padding: 0.5rem 0.75rem;
      position: fixed;
      right: 1rem;
      z-index: 1000;
    }}
  </style>
</head>
<body>
  <iframe
    id="seq-lm-report-frame"
    src={json.dumps(snapshot_name)}
    title="seq_lm analysis report"
  ></iframe>
  <div id="seq-lm-live-status" role="status"></div>
  <script>
    (() => {{
      const frame = document.getElementById("seq-lm-report-frame");
      const status = document.getElementById("seq-lm-live-status");
      const stateUrl = new URL({json.dumps(state_name)}, window.location.href);
      let displayedRevision = {json.dumps(report_revision)};
      let frameReady = false;
      let pendingViewState = null;
      let pendingRevision = null;

      const text = (element) =>
        element ? element.textContent.trim().replace(/\\s+/g, " ") : null;

      const tabKey = (tabList, index) =>
        tabList.dataset.seqLmTabKey || `tab-list-${{index}}`;

      const captureViewState = () => {{
        try {{
          const doc = frame.contentDocument;
          const win = frame.contentWindow;
          const tabs = Object.fromEntries(
            [...doc.querySelectorAll('[role="tablist"]')].map(
              (tabList, index) => [
                tabKey(tabList, index),
                text(tabList.querySelector('[data-bs-toggle="tab"].active'))
              ]
            )
          );
          return {{ tabs, scrollX: win.scrollX, scrollY: win.scrollY }};
        }} catch (error) {{
          return {{ tabs: {{}}, scrollX: 0, scrollY: 0 }};
        }}
      }};

      const restoreViewState = (viewState) => {{
        if (!viewState) return;
        try {{
          const doc = frame.contentDocument;
          const win = frame.contentWindow;
          const tabLists = [...doc.querySelectorAll('[role="tablist"]')];
          const tabListsByKey = new Map(
            tabLists.map((tabList, index) => [tabKey(tabList, index), tabList])
          );
          Object.entries(viewState.tabs || {{}}).forEach(([key, label]) => {{
            const tabList = tabListsByKey.get(key);
            if (!label || !tabList) return;
            const button = [...tabList.querySelectorAll(
              '[data-bs-toggle="tab"]'
            )].find((candidate) => text(candidate) === label);
            if (button) button.click();
          }});
          win.setTimeout(
            () => win.dispatchEvent(new Event("resize")),
            0
          );
          window.setTimeout(
            () => win.scrollTo(viewState.scrollX, viewState.scrollY),
            100
          );
        }} catch (error) {{
          console.warn("Unable to restore live report view state.", error);
        }}
      }};

      const showStatus = (message) => {{
        status.textContent = message;
        status.style.display = message ? "block" : "none";
      }};

      frame.addEventListener("load", () => {{
        try {{
          frameReady = Boolean(
            frame.contentDocument && frame.contentDocument.querySelector("main")
          );
        }} catch (error) {{
          frameReady = false;
        }}
        if (pendingRevision !== null) {{
          restoreViewState(pendingViewState);
          displayedRevision = pendingRevision;
          pendingViewState = null;
          pendingRevision = null;
        }}
        showStatus(
          frameReady ? "" : "Waiting for update."
        );
      }});

      const refreshReport = async () => {{
        try {{
          const requestUrl = new URL(stateUrl);
          requestUrl.searchParams.set("_", Date.now().toString());
          const response = await fetch(requestUrl, {{ cache: "no-store" }});
          if (!response.ok) {{
            throw new Error(`state request returned ${{response.status}}`);
          }}
          const state = await response.json();
          const nextRevision = state.report_revision;
          if (!Number.isSafeInteger(nextRevision) || nextRevision < 0) {{
            throw new Error("invalid report revision");
          }}
          const displayedNumber = Number(displayedRevision);
          const nextNumber = Number(nextRevision);
          const stateIsOlder = (
            Number.isFinite(displayedNumber) &&
            Number.isFinite(nextNumber) &&
            nextNumber < displayedNumber
          );
          if (
            stateIsOlder ||
            pendingRevision !== null ||
            (nextRevision === displayedRevision && frameReady)
          ) {{
            showStatus(
              frameReady ? "" : "Waiting for update."
            );
            return;
          }}

          const snapshotUrl = new URL(state.snapshot, window.location.href);
          snapshotUrl.searchParams.set("revision", String(nextRevision));
          const snapshotResponse = await fetch(snapshotUrl, {{
            method: "HEAD",
            cache: "no-store"
          }});
          if (!snapshotResponse.ok) {{
            throw new Error(
              `snapshot request returned ${{snapshotResponse.status}}`
            );
          }}
          const expectedBytes = Number(state.snapshot_bytes);
          const contentLength = snapshotResponse.headers.get("content-length");
          const availableBytes = contentLength === null
            ? null
            : Number(contentLength);
          if (
            Number.isFinite(expectedBytes) &&
            availableBytes !== null &&
            Number.isFinite(availableBytes) &&
            expectedBytes !== availableBytes
          ) {{
            throw new Error("snapshot is still being published");
          }}

          pendingViewState = captureViewState();
          pendingRevision = nextRevision;
          frame.src = snapshotUrl.href;
        }} catch (error) {{
          showStatus("Waiting for update.");
          console.warn("Live report refresh deferred.", error);
        }}
      }};

      window.setInterval(refreshReport, {refresh_milliseconds});
    }})();
  </script>
</body>
</html>
"""


def write_report(
    report,
    report_path,
    latest_batch,
    refresh_seconds,
    subtitle_notice=None,
    report_revision=None,
):
    """Write either a static report or a live shell with a versioned snapshot."""
    if refresh_seconds <= 0:
        report.write(report_path)
        apply_report_branding(report_path, subtitle_notice)
        return

    if report_revision is None:
        report_revision = _next_report_revision(report_path)
    if isinstance(report_revision, bool) or not isinstance(report_revision, int) or report_revision < 0:
        raise ValueError("Report revision must be a nonnegative integer")
    if report_revision < _next_report_revision(report_path):
        raise ValueError("Refusing to overwrite or publish an older report revision")
    report_path, state_path, snapshot_path = _live_report_paths(report_path, report_revision)
    report.write(snapshot_path)
    apply_report_branding(snapshot_path, subtitle_notice)
    report_path.write_text(
        _live_report_shell(
            state_path.name,
            snapshot_path.name,
            report_revision,
            refresh_seconds,
        )
    )
    state_path.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "report_revision": report_revision,
                "latest_batch": str(latest_batch),
                "snapshot": snapshot_path.name,
                "snapshot_bytes": snapshot_path.stat().st_size,
            },
            separators=(",", ":"),
        )
        + "\n"
    )


def main(args):
    """Run the entry point."""
    if args.lfc_cutoff < 0:
        raise ValueError("--lfc-cutoff must be nonnegative.")
    if not 0 < args.padj_cutoff <= 1:
        raise ValueError("--padj-cutoff must be greater than 0 and at most 1.")
    if args.gene_set_enrichment and args.differential_results is None:
        raise ValueError(
            "--gene-set-enrichment requires --differential-results."
        )
    if args.temporal_analysis and not args.gene_set_enrichment:
        raise ValueError(
            "--temporal-analysis requires --gene-set-enrichment."
        )

    logger = get_named_logger("SeqLMReport")
    samples = load_qc_samples(args.samples)
    differential = None
    stability_results = None
    fry_results = None
    gsva_results = None
    temporal_results = None
    imodulon_results = None
    transcript_biotypes = args.transcript_biotypes
    if args.differential_results is not None:
        differential = load_differential_results(
            args.differential_results,
            samples.samples_df,
        )
        stability_results = load_stability_results(
            args.stability_results,
            differential.sample_metadata,
            args.stability_behavior,
        )
        if args.gene_set_enrichment:
            fry_results = load_gene_set_results(
                args.differential_results,
                differential,
            )
            gsva_results = load_gsva_results(
                args.differential_results,
                differential,
            )
            if args.temporal_analysis:
                temporal_results = load_temporal_results(
                    args.differential_results,
                    samples.samples_df,
                    differential,
                    gsva_results,
                )
    if args.imodulon_results is not None:
        imodulon_results = load_imodulon_results(
            args.imodulon_results,
            args.imodulon_batch,
            args.imodulon_sequence,
            args.imodulon_analysis_index,
            samples.samples_df,
        )

    report = labs_report(
        labs,
        "seq_lm analysis report",
        "qc_report",
        args.params,
        args.versions,
        "workflow",
    )

    with report.add_section("Analysis", "Analysis"):  # type: ignore
        primary_tabs = Tabs()
        with primary_tabs.add_tab("Quality Control"):
            qc_tabs = Tabs()
            with qc_tabs.add_tab("Read flow"):
                add_sample_read_fate_sankeys(samples.sample_results)
            if transcript_biotypes is not None:
                with qc_tabs.add_tab("Transcript biotypes"):
                    add_transcript_biotype_composition(
                        transcript_biotypes,
                        samples.sample_results,
                    )
            with qc_tabs.add_tab("Metrics"):
                DataTable.from_pandas(
                    create_nanoplot_metrics_table(samples.sample_results),
                    use_index=False,
                )
            with qc_tabs.add_tab("Read length"):
                add_sample_hists(
                    samples.sample_results,
                    x_column="lengths",
                    title="Read length histogram",
                    x_axis_label="Read length",
                    y_axis_label="Number of reads",
                )
            with qc_tabs.add_tab("Read quality"):
                add_sample_2d_kdes(
                    samples.sample_results,
                    x_column="lengths",
                    y_column="quals",
                    title="Read length vs Read quality",
                    x_axis_label="Read length",
                    y_axis_label="Read quality",
                )
            with qc_tabs.add_tab("Mapping quality"):
                add_sample_2d_kdes(
                    samples.sample_results,
                    x_column="lengths",
                    y_column="mapQ",
                    title="Read length vs Mapping quality",
                    x_axis_label="Read length",
                    y_axis_label="Mapping quality",
                )
            with qc_tabs.add_tab("Samples"):
                DataTable.from_pandas(samples.samples_df)

        if differential is not None:
            with primary_tabs.add_tab("Differential Analysis"):
                add_differential_analysis(
                    differential,
                    args.lfc_cutoff,
                    args.padj_cutoff,
                    stability_results,
                )

            if args.gene_set_enrichment:
                with primary_tabs.add_tab("Gene Set Enrichment"):
                    analysis_tabs = Tabs()
                    with analysis_tabs.add_tab("GSVA scores"):
                        add_gsva_scores(
                            gsva_results,
                            differential.condition_colors,
                        )
                    with analysis_tabs.add_tab("GSVA differential"):
                        add_gsva_differential(
                            gsva_results,
                            differential.condition_colors,
                            args.padj_cutoff,
                        )
                    with analysis_tabs.add_tab("fry enrichment"):
                        add_gene_set_enrichment(
                            fry_results,
                            differential.condition_colors,
                            args.padj_cutoff,
                        )

                if args.temporal_analysis:
                    with primary_tabs.add_tab("Temporal Analysis"):
                        add_temporal_analysis(
                            temporal_results,
                            differential.condition_colors,
                        )

        if imodulon_results is not None:
            with primary_tabs.add_tab("iModulon Analysis"):
                add_imodulon_analysis(
                    imodulon_results,
                )

    write_report(  # type: ignore
        report,
        args.report,
        args.latest_batch,
        args.refresh_seconds,
        args.dea_readiness_notice,
        args.report_revision,
    )
    logger.info(f"Analysis report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("qc_report")
    parser.add_argument("report", help="Analysis report output HTML file")
    parser.add_argument(
        "--samples",
        required=True,
        help="TSV containing the current QC report sample rows.",
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="Directory containing CSVs containing name,version.",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="JSON file containing the workflow parameter key/values.",
    )
    parser.add_argument(
        "--latest-batch",
        default="unknown",
        help="Latest live analysis batch index represented in the report.",
    )
    parser.add_argument("--report-revision", type=int, help="Persistent publication revision, independent of local batch indices.")
    parser.add_argument("--imodulon-results", help="Immutable ICA snapshot directory.")
    parser.add_argument("--imodulon-batch", type=int, help="Expected ICA source batch index.")
    parser.add_argument("--imodulon-sequence", type=int, help="Expected ICA report sequence.")
    parser.add_argument("--imodulon-analysis-index", type=int, help="Expected ICA analysis index.")
    parser.add_argument(
        "--differential-results",
        help=(
            "Differential-expression batch directory containing edgeR, fry, "
            "and GSVA outputs."
        ),
    )
    parser.add_argument(
        "--transcript-biotypes",
        help=(
            "TSV containing Oarfish estimated-read fractions for canonical "
            "transcript biotypes."
        ),
    )
    parser.add_argument(
        "--stability-results",
        help="Per-sample stability audit TSV for the represented DE snapshot.",
    )
    parser.add_argument(
        "--stability-behavior",
        choices=("disabled", "log", "terminate"),
        default="disabled",
        help="Configured DE stability behavior [default: %(default)s].",
    )
    parser.add_argument(
        "--gene-set-enrichment",
        action="store_true",
        help="Include fry and GSVA gene-set results in the report.",
    )
    parser.add_argument(
        "--temporal-analysis",
        action="store_true",
        help=(
            "Include descriptive elapsed-minute gene-set trajectories; requires "
            "--gene-set-enrichment and sample order metadata."
        ),
    )
    parser.add_argument(
        "--dea-readiness-notice",
        help="Show this DEA precondition notice and omit unavailable result tabs.",
    )
    parser.add_argument(
        "--lfc-cutoff",
        default=1.0,
        type=float,
        help="Absolute log2 fold-change cutoff for differential plots.",
    )
    parser.add_argument(
        "--padj-cutoff",
        default=0.05,
        type=float,
        help="edgeR FDR cutoff for differential plots.",
    )
    parser.add_argument(
        "--refresh-seconds",
        default=5,
        type=int,
        help="Browser auto-refresh interval; use 0 to disable.",
    )
    return parser
