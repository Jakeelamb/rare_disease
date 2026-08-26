"""Local-only extraction of the restricted phenotype table into typed HPO evidence."""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from docx import Document
from pydantic import BaseModel

from .provenance import sha256

HPO_PATTERN = re.compile(r"HP:\d{7}")


class PhenotypeObservation(BaseModel):
    feature: str
    hpo_term: str
    hpo_id: str
    notes: str
    subject: Literal["proband", "family", "unresolved"]
    present: bool = True


class PhenotypeManifest(BaseModel):
    schema_version: str = "1.1"
    source_name: str
    source_sha256: str
    observations: tuple[PhenotypeObservation, ...]
    curation_status: Literal["unreviewed", "reviewed"] = "unreviewed"
    reviewer: str | None = None
    reviewed_at: str | None = None


def normalize_hpo_id(value: str) -> str:
    compact = re.sub(r"\s+", "", value.upper())
    match = HPO_PATTERN.search(compact)
    if not match:
        raise ValueError(f"no valid HPO identifier in {value!r}")
    return match.group(0)


def extract_phenotype(document_path: Path) -> PhenotypeManifest:
    document = Document(str(document_path))
    target = None
    required = {"Clinical Feature", "HPO Term", "HPO ID", "Presentation / Notes"}
    for table in document.tables:
        headers = {cell.text.strip() for cell in table.rows[0].cells}
        if required <= headers:
            target = table
            break
    if target is None:
        raise ValueError("clinical phenotype table with required headers was not found")

    header_map = {cell.text.strip(): index for index, cell in enumerate(target.rows[0].cells)}
    observations: list[PhenotypeObservation] = []
    for row in target.rows[1:]:
        values = [cell.text.strip() for cell in row.cells]
        feature = values[header_map["Clinical Feature"]]
        notes = values[header_map["Presentation / Notes"]]
        observations.append(
            PhenotypeObservation(
                feature=feature,
                hpo_term=values[header_map["HPO Term"]],
                hpo_id=normalize_hpo_id(values[header_map["HPO ID"]]),
                notes=notes,
                subject="unresolved",
            )
        )
    if not observations:
        raise ValueError("phenotype table contains no observations")
    return PhenotypeManifest(
        source_name=document_path.name,
        source_sha256=sha256(document_path),
        observations=tuple(observations),
    )


def write_phenotype(manifest: PhenotypeManifest, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(manifest.model_dump_json(indent=2) + "\n", encoding="utf-8")


CURATION_DECISIONS = {
    "proband_present": ("proband", True),
    "proband_absent": ("proband", False),
    "family_present": ("family", True),
    "family_absent": ("family", False),
}


def curate_phenotype(
    raw: PhenotypeManifest, decisions: dict[int, str], *, reviewer: str
) -> PhenotypeManifest:
    expected = set(range(1, len(raw.observations) + 1))
    if set(decisions) != expected:
        raise ValueError("curation requires exactly one decision for every phenotype row")
    if not reviewer.strip():
        raise ValueError("phenotype curator name cannot be empty")
    observations: list[PhenotypeObservation] = []
    for index, observation in enumerate(raw.observations, start=1):
        decision = decisions[index]
        if decision not in CURATION_DECISIONS:
            raise ValueError(f"unsupported phenotype decision for row {index}: {decision}")
        subject, present = CURATION_DECISIONS[decision]
        observations.append(observation.model_copy(update={"subject": subject, "present": present}))
    curated = raw.model_copy(
        update={
            "observations": tuple(observations),
            "curation_status": "reviewed",
            "reviewer": reviewer.strip(),
            "reviewed_at": datetime.now(UTC).isoformat(),
        }
    )
    validate_phenotype_curation(curated)
    return curated


def validate_phenotype_curation(manifest: PhenotypeManifest) -> None:
    if manifest.curation_status != "reviewed" or not manifest.reviewer or not manifest.reviewed_at:
        raise ValueError("manual phenotype curation is required before gene scoring")
    if any(observation.subject == "unresolved" for observation in manifest.observations):
        raise ValueError("manual phenotype curation left unresolved rows")


def render_phenotype_review(manifest: PhenotypeManifest) -> str:
    lines = [
        "# Private phenotype extraction review",
        "",
        "> Assign every row explicitly before phenotype-to-gene scoring.",
        "",
        "| Row | Feature | HPO | Raw notes | Decision |",
        "|---:|---|---|---|---|",
    ]
    for index, observation in enumerate(manifest.observations, start=1):
        values = [observation.feature, observation.hpo_id, observation.notes]
        escaped = [value.replace("|", "\\|").replace("\n", " ") for value in values]
        decision = (
            "unresolved"
            if observation.subject == "unresolved"
            else f"{observation.subject}_{'present' if observation.present else 'absent'}"
        )
        lines.append(f"| {index} | {escaped[0]} | `{escaped[1]}` | {escaped[2]} | {decision} |")
    lines.extend(
        [
            "",
            "Allowed decisions: `proband_present`, `proband_absent`, "
            "`family_present`, `family_absent`.",
        ]
    )
    return "\n".join(lines) + "\n"


def public_summary(manifest: PhenotypeManifest) -> str:
    counts: dict[str, int] = {}
    for observation in manifest.observations:
        counts[observation.subject] = counts.get(observation.subject, 0) + 1
    return json.dumps(
        {
            "observation_count": len(manifest.observations),
            "subject_counts": counts,
            "source_sha256": manifest.source_sha256,
        },
        indent=2,
    )
