import csv
from pathlib import Path

from ..models.sample import Sample


class SampleSheetError(ValueError):
    pass


def parse_is_live(value: str | None, alias: str | None = None) -> bool:
    normalized = (value or "").strip().lower()
    if normalized in {"", "true"}:
        return True
    if normalized == "false":
        return False
    sample_context = f" for sample '{alias}'" if alias else ""
    raise SampleSheetError(
        f"Invalid is_live value '{value}'{sample_context}. "
        "Expected true, false, or a blank value."
    )


def _parse_order(value: str | None, alias: str) -> int | None:
    normalized = (value or "").strip()
    if not normalized:
        return None
    try:
        return int(normalized)
    except ValueError as error:
        raise SampleSheetError(
            f"Invalid order value '{value}' for sample '{alias}'. Expected an integer."
        ) from error


def parse_sample_sheet(
    path: str | Path, *, require_unique_aliases: bool = False
) -> list[Sample]:
    sample_sheet = Path(path).expanduser()
    if not sample_sheet.is_file():
        raise SampleSheetError(f"Sample sheet does not exist: {sample_sheet}")

    try:
        with sample_sheet.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or [])
            required_fields = {"alias", "group", "bam_dir"}
            missing_fields = sorted(required_fields - fields)
            if missing_fields:
                raise SampleSheetError(
                    "Sample sheet is missing required fields: "
                    + ", ".join(missing_fields)
                )

            removed_fields = sorted({"id", "type"} & fields)
            if removed_fields:
                raise SampleSheetError(
                    "Sample sheet contains removed fields: "
                    + ", ".join(removed_fields)
                )

            samples: list[Sample] = []
            for row_number, row in enumerate(reader, start=2):
                if None in row:
                    raise SampleSheetError(
                        f"Unexpected number of cells in row {row_number}."
                    )

                alias = (row.get("alias") or "").strip()
                group = (row.get("group") or "").strip()
                bam_dir_value = (row.get("bam_dir") or "").strip()
                if not alias:
                    raise SampleSheetError(
                        f"Sample sheet contains an empty alias in row {row_number}."
                    )
                if not group:
                    raise SampleSheetError(
                        f"Sample sheet contains an empty group for alias '{alias}'."
                    )
                if not bam_dir_value:
                    raise SampleSheetError(
                        f"Sample sheet contains an empty bam_dir for sample '{group}/{alias}'."
                    )

                bam_dir = Path(bam_dir_value).expanduser()
                if not bam_dir.exists():
                    raise SampleSheetError(
                        f"BAM directory does not exist for sample '{group}/{alias}': {bam_dir}"
                    )
                if not bam_dir.is_dir():
                    raise SampleSheetError(
                        f"BAM path is not a directory for sample '{group}/{alias}': {bam_dir}"
                    )

                samples.append(
                    Sample(
                        alias=alias,
                        group=group,
                        bam_dir=bam_dir,
                        is_live=parse_is_live(row.get("is_live"), alias),
                        order=_parse_order(row.get("order"), alias),
                    )
                )
    except (OSError, csv.Error) as error:
        raise SampleSheetError(f"Could not parse sample sheet {sample_sheet}: {error}") from error

    if not samples:
        raise SampleSheetError("Sample sheet must contain at least one row.")

    identities = [(sample.group, sample.alias) for sample in samples]
    duplicate_identities = sorted(
        {identity for identity in identities if identities.count(identity) > 1}
    )
    if duplicate_identities:
        duplicates = ", ".join(f"{group}/{alias}" for group, alias in duplicate_identities)
        raise SampleSheetError(
            f"Sample aliases must be unique within each group. Duplicates: {duplicates}"
        )

    if require_unique_aliases:
        aliases = [sample.alias for sample in samples]
        duplicates = sorted({alias for alias in aliases if aliases.count(alias) > 1})
        if duplicates:
            raise SampleSheetError(
                "Sample aliases must be globally unique when starting acquisitions. "
                f"Duplicates: {', '.join(duplicates)}"
            )

    return samples
