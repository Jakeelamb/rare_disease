"""Validate and bundle manually reviewed Track 1 submission artifacts."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .submission import validate_submission

REQUIRED_MARKDOWN_FIELDS = (
    "**Entrant/team:**",
    "**Model number:**",
    "## Method abstract",
    "## Automation and downstream manual review",
    "## Detailed computational approach",
    "## Data and software sources",
    "## Generative-AI assistance and data handling",
    "**Required disclosure:**",
    "## Runtime, hardware, and cost",
    "## Secondary or incidental findings",
    "## Strengths, limitations, and interpretation boundary",
)


@dataclass(frozen=True)
class ReportValidation:
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    abstract_words: int | None

    @property
    def ok(self) -> bool:
        return not self.errors


@dataclass(frozen=True)
class SubmissionPackage:
    bundle_dir: Path
    manifest_path: Path
    checksums_path: Path
    manifest: dict[str, Any]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _method_abstract(text: str) -> str | None:
    match = re.search(r"^## Method abstract\s*$\n(.*?)(?=^##\s)", text, re.MULTILINE | re.DOTALL)
    return match.group(1).strip() if match else None


def validate_methods_report(path: Path) -> ReportValidation:
    errors: list[str] = []
    warnings: list[str] = []
    abstract_words: int | None = None

    if path.suffix.lower() not in {".md", ".pdf"}:
        errors.append("methods report must be Markdown (.md) or PDF (.pdf)")
        return ReportValidation(tuple(errors), tuple(warnings), abstract_words)
    if not path.is_file() or path.stat().st_size == 0:
        errors.append("methods report must be a non-empty file")
        return ReportValidation(tuple(errors), tuple(warnings), abstract_words)
    if path.suffix.lower() == ".pdf":
        warnings.append(
            "PDF contents were not semantically checked; use Markdown for full preflight"
        )
        return ReportValidation(tuple(errors), tuple(warnings), abstract_words)

    text = path.read_text(encoding="utf-8")
    for field in REQUIRED_MARKDOWN_FIELDS:
        if field not in text:
            errors.append(f"methods report is missing required field or section: {field}")

    abstract = _method_abstract(text)
    if abstract is None:
        errors.append("methods report must contain a bounded Method abstract section")
    else:
        abstract_words = len(re.findall(r"\b[\w'-]+\b", abstract))
        if abstract_words == 0:
            errors.append("method abstract must not be empty")
        elif abstract_words > 500:
            errors.append(f"method abstract contains {abstract_words} words; maximum is 500")

    if re.search(r"could not be independently\s+verified", text):
        warnings.append("AI account data-control setting is disclosed as unverified")
    if "phase is unresolved" not in text and "unresolved short-read phase" not in text:
        warnings.append("report does not use the expected unresolved-phase caution")

    return ReportValidation(tuple(errors), tuple(warnings), abstract_words)


def _artifact(path: Path, media_type: str) -> dict[str, str | int]:
    return {
        "filename": path.name,
        "media_type": media_type,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _write_immutable(path: Path, content: bytes) -> None:
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"refusing to overwrite different existing bundle file: {path}")
        return
    path.write_bytes(content)


def create_submission_package(
    *,
    csv_path: Path,
    report_path: Path,
    pixi_lock_path: Path,
    output_root: Path,
    repository_url: str,
    repository_commit: str,
    space_revision: str,
    dataset_revision: str,
) -> SubmissionPackage:
    csv_validation = validate_submission(csv_path)
    if not csv_validation.ok:
        raise ValueError("invalid submission CSV: " + "; ".join(csv_validation.errors))
    report_validation = validate_methods_report(report_path)
    if not report_validation.ok:
        raise ValueError("invalid methods report: " + "; ".join(report_validation.errors))
    if not pixi_lock_path.is_file():
        raise ValueError(f"Pixi lockfile does not exist: {pixi_lock_path}")
    if not re.fullmatch(r"[0-9a-f]{40}", repository_commit):
        raise ValueError("repository commit must be a full 40-character lowercase Git SHA")
    if not repository_url.startswith("https://github.com/"):
        raise ValueError("repository URL must start with https://github.com/")

    artifacts = [
        _artifact(csv_path, "text/csv"),
        _artifact(report_path, "text/markdown"),
        _artifact(pixi_lock_path, "application/toml"),
    ]
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "track": "Track 1 - Variant Prediction",
        "challenge_revisions": {
            "dataset": dataset_revision,
            "space": space_revision,
        },
        "repository": {
            "url": repository_url,
            "commit": repository_commit,
        },
        "environment": {
            "manager": "Pixi",
            "lockfile": pixi_lock_path.name,
            "lockfile_sha256": artifacts[2]["sha256"],
        },
        "submission_artifacts": artifacts[:2],
        "validation": {
            "csv": "pass",
            "csv_warnings": list(csv_validation.warnings),
            "report": "pass",
            "report_abstract_words": report_validation.abstract_words,
            "report_warnings": list(report_validation.warnings),
        },
    }
    identity = json.dumps(
        {
            "challenge_revisions": manifest["challenge_revisions"],
            "repository": manifest["repository"],
            "environment": manifest["environment"],
            "submission_artifacts": manifest["submission_artifacts"],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    bundle_id = hashlib.sha256(identity).hexdigest()[:12]
    bundle_dir = output_root / f"track1-{bundle_id}"
    bundle_dir.mkdir(parents=True, exist_ok=True)

    for source in (csv_path, report_path, pixi_lock_path):
        _write_immutable(bundle_dir / source.name, source.read_bytes())

    manifest_path = bundle_dir / "submission-manifest.json"
    manifest_content = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    _write_immutable(manifest_path, manifest_content)

    checksum_entries = [*artifacts, _artifact(manifest_path, "application/json")]
    checksums_path = bundle_dir / "SHA256SUMS"
    checksums_content = "".join(
        f"{item['sha256']}  {item['filename']}\n" for item in checksum_entries
    ).encode()
    _write_immutable(checksums_path, checksums_content)
    return SubmissionPackage(bundle_dir, manifest_path, checksums_path, manifest)
