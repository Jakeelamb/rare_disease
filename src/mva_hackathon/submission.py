"""Read-only checks for a manually prepared Track 1 CSV."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

SUBMISSION_COLUMNS = (
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
)


@dataclass(frozen=True)
class ValidationResult:
    errors: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.errors


def validate_submission(path: Path) -> ValidationResult:
    errors: list[str] = []
    warnings: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SUBMISSION_COLUMNS:
            errors.append(f"columns must exactly equal {','.join(SUBMISSION_COLUMNS)}")
        rows = list(reader)

    if not rows:
        errors.append("submission must contain at least one row")
    if len(rows) > 10:
        errors.append(f"submission contains {len(rows)} rows; maximum is 10")

    probabilities: list[float] = []
    for line_number, row in enumerate(rows, start=2):
        if row.get("proband_id") != "PROBAND01":
            errors.append(f"line {line_number}: proband_id must be PROBAND01")
        for suffix in ("1", "2"):
            chrom = (row.get(f"chrom_{suffix}") or "").strip()
            pos = (row.get(f"pos_{suffix}") or "").strip()
            ref = (row.get(f"ref_{suffix}") or "").strip()
            alt = (row.get(f"alt_{suffix}") or "").strip()
            present = [bool(chrom), bool(pos), bool(ref), bool(alt)]
            if suffix == "1" and not any(present):
                errors.append(f"line {line_number}: variant 1 is required")
            if any(present) and not all(present):
                errors.append(f"line {line_number}: variant {suffix} fields are incomplete")
            if chrom and not chrom.startswith("chr"):
                errors.append(f"line {line_number}: {chrom} must use the literal chr prefix")
            if pos:
                try:
                    if int(pos) < 1 or str(int(pos)) != pos:
                        raise ValueError
                except ValueError:
                    errors.append(f"line {line_number}: pos_{suffix} must be a positive integer")
            if ref != ref.upper() or alt != alt.upper():
                errors.append(f"line {line_number}: REF/ALT must be uppercase")
        try:
            probability = float(row.get("epcr") or "")
            if not 0 < probability <= 1:
                raise ValueError
            probabilities.append(probability)
        except ValueError:
            errors.append(f"line {line_number}: epcr must be numeric in (0,1]")
        if row.get("finding_type") not in {"primary", "secondary"}:
            errors.append(f"line {line_number}: finding_type must be primary or secondary")

    if probabilities != sorted(probabilities, reverse=True):
        errors.append("EPCR values must be sorted strictly descending")
    if len(set(probabilities)) != len(probabilities):
        errors.append("EPCR values must be unique to make evaluator ordering deterministic")
    if rows and not (rows[0].get("chrom_2") or "").strip():
        warnings.append("rank-1 hypothesis is not a compound-heterozygous pair")
    if any(row.get("finding_type") == "secondary" for row in rows):
        warnings.append("official scorer ignores finding_type; secondary rows still enter F-max")
    return ValidationResult(errors=tuple(errors), warnings=tuple(warnings))
