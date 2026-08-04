"""Compatibility helpers for ezCharts report generation."""
import inspect

from bokeh.resources import Resources


def patch_bokeh_resources():
    """Provide the Bokeh resource API and widget bundle needed by ezCharts."""
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
