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
)
from .qc_report_types.gene_set_plots import (
    add_gene_set_enrichment,
    load_gene_set_results,
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


def _live_report_paths(report_path, latest_batch):
    """Return stable-shell, state, and immutable snapshot paths."""
    report_path = Path(report_path)
    batch_token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(latest_batch)).strip("._")
    batch_token = batch_token or "unknown"
    snapshot_path = report_path.with_name(
        f"{report_path.stem}_snapshot_{batch_token}{report_path.suffix}"
    )
    state_path = report_path.with_name(f"{report_path.stem}_state.json")
    return report_path, state_path, snapshot_path


def _live_report_shell(
    state_name,
    snapshot_name,
    latest_batch,
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
      let displayedBatch = {json.dumps(str(latest_batch))};
      let frameReady = false;
      let pendingViewState = null;
      let pendingBatch = null;

      const text = (element) =>
        element ? element.textContent.trim().replace(/\\s+/g, " ") : null;

      const captureViewState = () => {{
        try {{
          const doc = frame.contentDocument;
          const win = frame.contentWindow;
          const tabs = [...doc.querySelectorAll('[role="tablist"]')].map(
            (tabList) => text(
              tabList.querySelector('[data-bs-toggle="tab"].active')
            )
          );
          return {{ tabs, scrollX: win.scrollX, scrollY: win.scrollY }};
        }} catch (error) {{
          return {{ tabs: [], scrollX: 0, scrollY: 0 }};
        }}
      }};

      const restoreViewState = (viewState) => {{
        if (!viewState) return;
        try {{
          const doc = frame.contentDocument;
          const win = frame.contentWindow;
          const tabLists = [...doc.querySelectorAll('[role="tablist"]')];
          viewState.tabs.forEach((label, index) => {{
            if (!label || !tabLists[index]) return;
            const button = [...tabLists[index].querySelectorAll(
              '[data-bs-toggle="tab"]'
            )].find((candidate) => text(candidate) === label);
            if (button) button.click();
          }});
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
        if (pendingBatch !== null) {{
          restoreViewState(pendingViewState);
          displayedBatch = pendingBatch;
          pendingViewState = null;
          pendingBatch = null;
        }}
        showStatus(
          frameReady ? "" : "Waiting for the next complete report update…"
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
          const nextBatch = String(state.latest_batch);
          const displayedNumber = Number(displayedBatch);
          const nextNumber = Number(nextBatch);
          const stateIsOlder = (
            Number.isFinite(displayedNumber) &&
            Number.isFinite(nextNumber) &&
            nextNumber < displayedNumber
          );
          if (
            stateIsOlder ||
            pendingBatch !== null ||
            (nextBatch === displayedBatch && frameReady)
          ) {{
            showStatus(
              frameReady ? "" : "Waiting for the next complete report update…"
            );
            return;
          }}

          const snapshotUrl = new URL(state.snapshot, window.location.href);
          snapshotUrl.searchParams.set("batch", nextBatch);
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
          pendingBatch = nextBatch;
          frame.src = snapshotUrl.href;
        }} catch (error) {{
          showStatus("Waiting for the next complete report update…");
          console.warn("Live report refresh deferred.", error);
        }}
      }};

      window.setInterval(refreshReport, {refresh_milliseconds});
    }})();
  </script>
</body>
</html>
"""


def write_report(report, report_path, latest_batch, refresh_seconds):
    """Write either a static report or a live shell with a versioned snapshot."""
    if refresh_seconds <= 0:
        report.write(report_path)
        apply_report_branding(report_path)
        return

    report_path, state_path, snapshot_path = _live_report_paths(
        report_path,
        latest_batch,
    )
    report.write(snapshot_path)
    apply_report_branding(snapshot_path)
    report_path.write_text(
        _live_report_shell(
            state_path.name,
            snapshot_path.name,
            latest_batch,
            refresh_seconds,
        )
    )
    state_path.write_text(
        json.dumps(
            {
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

    logger = get_named_logger("SeqLMReport")
    samples = load_qc_samples(args.samples)
    differential = load_differential_results(
        args.differential_results,
        samples.samples_df,
    )
    gene_sets = load_gene_set_results(
        args.differential_results,
        differential,
    )

    report = labs_report(
        labs,
        "seq_lm analysis report",
        "qc_report",
        args.params,
        args.versions,
        "workflow",
    )

    with report.add_section("Quality Control", "Quality Control"):  # type: ignore
        tabs = Tabs()
        with tabs.add_tab("Read flow"):
            add_sample_read_fate_sankeys(samples.sample_results)
        with tabs.add_tab("Metrics"):
            DataTable.from_pandas(
                create_nanoplot_metrics_table(samples.sample_results),
                use_index=False,
            )
        with tabs.add_tab("Read length"):
            add_sample_hists(
                samples.sample_results,
                x_column="lengths",
                title="Read length histogram",
                x_axis_label="Read length",
                y_axis_label="Number of reads",
            )
        with tabs.add_tab("Read quality"):
            add_sample_2d_kdes(
                samples.sample_results,
                x_column="lengths",
                y_column="quals",
                title="Read length vs Read quality",
                x_axis_label="Read length",
                y_axis_label="Read quality",
            )
        with tabs.add_tab("Mapping quality"):
            add_sample_2d_kdes(
                samples.sample_results,
                x_column="lengths",
                y_column="mapQ",
                title="Read length vs Mapping quality",
                x_axis_label="Read length",
                y_axis_label="Mapping quality",
            )
        with tabs.add_tab("Samples"):
            DataTable.from_pandas(samples.samples_df)

    with report.add_section(  # type: ignore
        "Differential Analysis",
        "Differential Analysis",
    ):
        add_differential_analysis(
            differential,
            args.lfc_cutoff,
            args.padj_cutoff,
        )

    with report.add_section(  # type: ignore
        "Gene Set Enrichment",
        "Gene Set Enrichment",
    ):
        add_gene_set_enrichment(
            gene_sets,
            differential.condition_colors,
            args.padj_cutoff,
        )

    write_report(  # type: ignore
        report,
        args.report,
        args.latest_batch,
        args.refresh_seconds,
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
    parser.add_argument(
        "--differential-results",
        required=True,
        help=(
            "edgeR batch result directory containing feature counts, "
            "contrast results, and fry outputs."
        ),
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
