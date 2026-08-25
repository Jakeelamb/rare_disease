import csv
from pathlib import Path

from mva_hackathon.submission import SUBMISSION_COLUMNS, validate_submission


def test_validator_rejects_literal_contract_hazards(tmp_path: Path) -> None:
    output = tmp_path / "bad.csv"
    output.write_text(
        ",".join(SUBMISSION_COLUMNS)
        + "\nPROBAND01,15,100,A,G,,,,,0.9,primary,missing chr prefix\n",
        encoding="utf-8",
    )

    result = validate_submission(output)

    assert not result.ok
    assert any("literal chr prefix" in error for error in result.errors)
    assert any("rank-1 hypothesis" in warning for warning in result.warnings)


def test_validator_rejects_tied_epcr(tmp_path: Path) -> None:
    output = tmp_path / "tied.csv"
    rows = [
        ["PROBAND01", "chr1", "100", "A", "G", "", "", "", "", "0.5", "primary", ""],
        ["PROBAND01", "chr2", "200", "C", "T", "", "", "", "", "0.5", "primary", ""],
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(SUBMISSION_COLUMNS)
        writer.writerows(rows)

    result = validate_submission(output)

    assert any("unique" in error for error in result.errors)


def test_validator_requires_complete_primary_variant_and_integer_positions(tmp_path: Path) -> None:
    output = tmp_path / "invalid-primary.csv"
    rows = [
        ["PROBAND01", "", "", "", "", "", "", "", "", "0.9", "primary", ""],
        [
            "PROBAND01",
            "chr1",
            "1.5",
            "A",
            "G",
            "",
            "",
            "",
            "",
            "0.8",
            "primary",
            "",
        ],
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(SUBMISSION_COLUMNS)
        writer.writerows(rows)

    result = validate_submission(output)

    assert any("variant 1 is required" in error for error in result.errors)
    assert any("pos_1 must be a positive integer" in error for error in result.errors)
