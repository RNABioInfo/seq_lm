"""Compatibility helpers for ezCharts report generation."""
import inspect

from bokeh.resources import Resources


def patch_bokeh_resources():
    """Provide the Bokeh resource API expected by ezCharts releases in use."""
    if not hasattr(Resources, "components_for") and hasattr(Resources, "components"):
        Resources.components_for = Resources.components


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
