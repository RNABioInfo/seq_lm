"""Contrasts select analyses independently of the full QC cohort."""

import csv
from pathlib import Path

import pytest

from workflow_glue.gene_contrasts import read_gene_contrasts, resolve_ica_groups
from workflow_glue import check_sample_sheet


def test_order_quoting_and_reversed_contrasts(tmp_path):
    path = tmp_path / "contrasts.csv"
    path.write_text('treated,baseline\n"quoted, group", treated\nbaseline,treated\n')
    assert read_gene_contrasts(path, {"baseline", "treated", "quoted, group", "QC only"}) == [
        ("treated", "baseline"), ("quoted, group", "treated"), ("baseline", "treated")]


def test_quoted_fields_allow_surrounding_spaces_and_literal_quotes(tmp_path):
    path = tmp_path / "contrasts.csv"
    path.write_text('  "quoted, "" group"  ,  "baseline"  \n')
    assert read_gene_contrasts(path, {'quoted, " group', 'baseline'}) == [('quoted, " group', 'baseline')]


@pytest.mark.parametrize("contents", ["", "\n", "A\n", "A,B,C\n", "A,\n", "A,A\n", "A,X\n", "A,B\nA,B\n", '"A,B\n'])
def test_rejects_invalid_csv(tmp_path, contents):
    path = tmp_path / "contrasts.csv"
    path.write_text(contents)
    with pytest.raises((ValueError, csv.Error)):
        read_gene_contrasts(path, {"A", "B"})


def test_checker_accepts_no_control_with_contrasts_and_validates_ica_independently(tmp_path):
    sheet = tmp_path / "samples.csv"
    sheet.write_text("alias,group,bam_dir\n" + "".join(
        f"{group}_{i},{group},{tmp_path}\n" for group in ["A", "B", "QC"] for i in [1, 2]))
    contrasts = tmp_path / "contrasts.csv"
    contrasts.write_text("A,B\n")
    args = check_sample_sheet.argparser().parse_args([
        str(sheet), "--gene_contrast_file", str(contrasts), "--ica_analysis",
        "--ica_control_group", "B", "--ica_condition_groups", "A"])
    check_sample_sheet.main(args)
    args.ica_control_group = "missing"
    with pytest.raises(SystemExit):
        check_sample_sheet.main(args)
    args.gene_contrast_file = None
    with pytest.raises(SystemExit):
        check_sample_sheet.main(args)
    args.differential_expression = False
    args.ica_control_group = "B"
    check_sample_sheet.main(args)


def test_ica_selection_includes_control_and_defaults_to_all():
    groups = {"reference", "treated", "QC"}
    assert resolve_ica_groups(groups, "reference", " treated ") == {"reference", "treated"}
    assert resolve_ica_groups(groups, "reference", None) == groups
    assert resolve_ica_groups(groups, "reference", " ") == groups
    for invalid in ["missing", "treated,", "treated,treated"]:
        with pytest.raises(ValueError):
            resolve_ica_groups(groups, "reference", invalid)
