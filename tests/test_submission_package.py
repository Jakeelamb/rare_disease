import json
from pathlib import Path

import pytest

from mva_hackathon.submission_package import (
    create_submission_package,
    sha256_file,
    validate_methods_report,
)


def _report(abstract: str = "A transparent method with explicit limitations.") -> str:
    return f"""# Track 1 variant prediction and methods report

**Entrant/team:** Test Team
**Model number:** 1

## Method abstract

{abstract}

## Automation and downstream manual review

Reviewed manually.

## Detailed computational approach

Detailed method.

## Data and software sources

Public sources.

## Generative-AI assistance and data handling

**Required disclosure:** No AI used.

## Runtime, hardware, and cost

One minute, zero marginal cost.

## Secondary or incidental findings

None.

## Strengths, limitations, and interpretation boundary

Trans phase is unresolved.
"""


def _csv() -> str:
    return (
        "proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,"
        "epcr,finding_type,notes\n"
        "PROBAND01,chr15,100,A,G,chr15,200,C,T,0.9,primary,test pair\n"
    )


def test_report_preflight_checks_template_sections_and_abstract_limit(tmp_path: Path) -> None:
    report = tmp_path / "report.md"
    report.write_text(_report(), encoding="utf-8")

    result = validate_methods_report(report)

    assert result.ok
    assert result.abstract_words == 6

    report.write_text(_report("word " * 501), encoding="utf-8")
    result = validate_methods_report(report)
    assert not result.ok
    assert any("maximum is 500" in error for error in result.errors)


def test_submission_package_is_content_addressed_and_idempotent(tmp_path: Path) -> None:
    csv_path = tmp_path / "entry.csv"
    report_path = tmp_path / "report.md"
    lock_path = tmp_path / "pixi.lock"
    csv_path.write_text(_csv(), encoding="utf-8")
    report_path.write_text(_report(), encoding="utf-8")
    lock_path.write_text("version = 6\n", encoding="utf-8")

    kwargs = {
        "csv_path": csv_path,
        "report_path": report_path,
        "pixi_lock_path": lock_path,
        "output_root": tmp_path / "packages",
        "repository_url": "https://github.com/example/repo",
        "repository_commit": "a" * 40,
        "space_revision": "b" * 40,
        "dataset_revision": "c" * 40,
    }
    first = create_submission_package(**kwargs)
    second = create_submission_package(**kwargs)

    assert first.bundle_dir == second.bundle_dir
    assert first.bundle_dir.name.startswith("track1-")
    assert (first.bundle_dir / csv_path.name).read_bytes() == csv_path.read_bytes()
    assert (first.bundle_dir / report_path.name).read_bytes() == report_path.read_bytes()
    assert (first.bundle_dir / lock_path.name).read_bytes() == lock_path.read_bytes()
    manifest = json.loads(first.manifest_path.read_text(encoding="utf-8"))
    assert manifest["repository"]["commit"] == "a" * 40
    assert manifest["validation"]["csv"] == "pass"
    assert manifest["validation"]["report"] == "pass"
    assert sha256_file(first.manifest_path) in first.checksums_path.read_text(encoding="utf-8")


def test_submission_package_rejects_invalid_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "entry.csv"
    report_path = tmp_path / "report.md"
    lock_path = tmp_path / "pixi.lock"
    csv_path.write_text("bad,header\n", encoding="utf-8")
    report_path.write_text(_report(), encoding="utf-8")
    lock_path.write_text("version = 6\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid submission CSV"):
        create_submission_package(
            csv_path=csv_path,
            report_path=report_path,
            pixi_lock_path=lock_path,
            output_root=tmp_path / "packages",
            repository_url="https://github.com/example/repo",
            repository_commit="a" * 40,
            space_revision="b" * 40,
            dataset_revision="c" * 40,
        )
