"""Compatibility helpers for ezCharts report generation."""
import inspect
from pathlib import Path
import re

from bokeh.resources import Resources


_BRAND_COLOR = "#004191"
# Keep the report portable and container-safe by embedding the supplied PNG.
_BRAND_LOGO_DATA_URI = (
    "data:image/png;base64,"
    "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"
    "JahJCu5V+KySAPCbiBWwAAAABJRU5ErkJggg=="
)
_EPI2ME_HEADER_LINK = re.compile(
    r'<a\b[^>]*\bhref=(?P<quote>["\'])https://labs\.epi2me\.io/?'
    r'(?P=quote)[^>]*>.*?</a>',
    flags=re.IGNORECASE | re.DOTALL,
)
_BRAND_STYLES = f"""
    <style id="seq-lm-report-branding">
      header .bg-dark {{
        background-color: {_BRAND_COLOR} !important;
      }}
      .seq-lm-report-brand {{
        background: #fff;
        border-radius: 0.375rem;
        padding: 0.35rem 0.6rem;
      }}
      .seq-lm-report-brand img {{
        display: block;
        height: 38px;
        max-width: min(46vw, 180px);
        object-fit: contain;
        width: auto;
      }}
      main .tab-content.p-5 {{
        padding: 1.25rem !important;
      }}
      main .tab-content .tab-content.p-5 {{
        padding: 0.75rem 0.25rem !important;
      }}
      main .tab-content .tab-content .tab-content.p-5 {{
        padding: 0.5rem 0 !important;
      }}
      .seq-lm-primary-tablist {{
        border-bottom: 2px solid {_BRAND_COLOR};
        gap: 0.35rem;
      }}
      .seq-lm-primary-tablist > .nav-item > .nav-link {{
        border: 1px solid transparent;
        border-radius: 0.375rem 0.375rem 0 0;
        color: #495057 !important;
        font-weight: 600;
        margin-right: 0 !important;
        padding: 0.75rem 1rem !important;
      }}
      .seq-lm-primary-tablist > .nav-item > .nav-link.active {{
        background: {_BRAND_COLOR};
        border-color: {_BRAND_COLOR};
        color: #fff !important;
      }}
      .seq-lm-subtablist {{
        border-bottom-color: #adb5bd;
        margin-top: 0.35rem;
      }}
      .seq-lm-subtablist > .nav-item > .nav-link {{
        padding-bottom: 0.4rem !important;
        padding-top: 0.4rem !important;
      }}
      @media (max-width: 767.98px) {{
        .seq-lm-primary-tablist > .nav-item > .nav-link {{
          padding: 0.55rem 0.65rem !important;
        }}
        main .tab-content.p-5 {{
          padding: 0.75rem 0.25rem !important;
        }}
      }}
    </style>
"""
_BRAND_MARKUP = (
    '<div class="seq-lm-report-brand d-flex align-items-center '
    'mb-md-0 me-md-auto" aria-label="RNA BioInfo AUCG">'
    f'<img src="{_BRAND_LOGO_DATA_URI}" alt="RNA BioInfo AUCG logo">'
    "</div>"
)
_REPORT_NAVIGATION_SCRIPT = r"""
    <script id="seq-lm-report-navigation">
      (() => {
        const text = (element) => element
          ? element.textContent.trim().replace(/\s+/g, " ")
          : "";

        const directTabButtons = (tabList) => [...tabList.children]
          .flatMap((child) => [...child.querySelectorAll(":scope > button")])
          .filter((button) => button.matches('[data-bs-toggle="tab"]'));

        const paneLabel = (pane) => {
          const labelledBy = (pane.getAttribute("aria-labelledby") || "")
            .replace(/^#/, "");
          return text(document.getElementById(labelledBy));
        };

        const markTabLevels = () => {
          const tabLists = [...document.querySelectorAll('[role="tablist"]')];
          tabLists.forEach((tabList) => {
            const labels = directTabButtons(tabList).map(text);
            const isPrimary = [
              "Quality Control",
              "Differential Analysis",
              "Gene Set Enrichment"
            ].every((label) => labels.includes(label));
            tabList.classList.toggle("seq-lm-primary-tablist", isPrimary);
            tabList.classList.toggle("seq-lm-subtablist", !isPrimary);

            const ancestors = [];
            let parent = tabList.parentElement;
            while (parent) {
              if (parent.classList && parent.classList.contains("tab-pane")) {
                const label = paneLabel(parent);
                if (label) ancestors.unshift(label);
              }
              parent = parent.parentElement;
            }
            const ownIdentity = labels.length
              ? labels.join(" / ")
              : text(tabList.querySelector(".dropdown-toggle"));
            tabList.dataset.seqLmTabKey = [...ancestors, ownIdentity].join(" > ");
          });
        };

        const resizeVisibleCharts = () => {
          window.dispatchEvent(new Event("resize"));
          if (window.Bokeh && window.Bokeh.index) {
            Object.values(window.Bokeh.index).forEach((view) => {
              if (view && typeof view.resize_layout === "function") {
                view.resize_layout();
              }
            });
          }
          if (window.echarts) {
            document.querySelectorAll("[_echarts_instance_]").forEach((element) => {
              const chart = window.echarts.getInstanceByDom(element);
              if (chart) chart.resize();
            });
          }
        };

        markTabLevels();
        document.addEventListener("shown.bs.tab", () => {
          window.setTimeout(resizeVisibleCharts, 0);
          window.setTimeout(resizeVisibleCharts, 150);
        });
        window.addEventListener("load", () => {
          markTabLevels();
          window.setTimeout(resizeVisibleCharts, 0);
        });
      })();
    </script>
"""


def _insert_before_last(html, closing_tag, content):
    """Insert content before the document-level closing tag.

    Vendored JavaScript can contain strings such as ``</body>``. Using a
    first-match replacement would insert report markup into those strings.
    """
    index = html.rfind(closing_tag)
    if index == -1:
        return html
    return f"{html[:index]}{content}{html[index:]}"


def apply_report_branding(report_path):
    """Apply branding and robust nested-tab behavior to a LabsReport."""
    report_path = Path(report_path)
    html = report_path.read_text()
    html = _EPI2ME_HEADER_LINK.sub(_BRAND_MARKUP, html, count=1)
    if 'id="seq-lm-report-branding"' not in html and "</head>" in html:
        html = _insert_before_last(html, "</head>", _BRAND_STYLES)
    if 'id="seq-lm-report-navigation"' not in html and "</body>" in html:
        html = _insert_before_last(
            html,
            "</body>",
            _REPORT_NAVIGATION_SCRIPT,
        )
    report_path.write_text(html)


def patch_bokeh_resources():
    """Provide the chart resources needed by self-contained ezCharts reports."""
    if not hasattr(Resources, "components_for") and hasattr(Resources, "components"):
        Resources.components_for = Resources.components

    # ezCharts 0.16 only inlines the core BokehJS bundle. A report containing a
    # Bokeh widget (for example, the gene-set selector) then serializes a Select
    # model that BokehJS cannot resolve. Because ezCharts embeds every Bokeh plot
    # in one document, that one unresolved model prevents all Bokeh plots from
    # rendering, including the read QC and differential-analysis plots.
    from ezcharts.layout import resource as ezcharts_resource

    if getattr(ezcharts_resource.bokeh_js.func, "_seq_lm_widgets", False):
        return

    def get_bokeh_js_with_widgets():
        inline = ezcharts_resource.bk_inline
        components = inline.components_for("js")
        component_indices = {
            component: index for index, component in enumerate(components)
        }
        required = ("bokeh", "bokeh-widgets")
        if not all(component in component_indices for component in required):
            return ezcharts_resource.get_bokeh_js()
        javascript = "\n".join(
            inline.js_raw[component_indices[component]] for component in required
        )
        return ezcharts_resource.raw(javascript)

    get_bokeh_js_with_widgets._seq_lm_widgets = True
    ezcharts_resource.bokeh_js.func = get_bokeh_js_with_widgets


def labs_report(
    labs_module,
    title,
    workflow_name,
    params_path,
    versions_path,
    version,
):
    """Create a LabsReport across ezCharts constructor variants."""
    patch_bokeh_resources()

    args = [title, workflow_name, params_path, versions_path]
    if "workflow_version" in inspect.signature(labs_module.LabsReport).parameters:
        args.append(version)
    return labs_module.LabsReport(*args)
