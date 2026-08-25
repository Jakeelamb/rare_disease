"""Input identity and environment provenance."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def tool_version(command: list[str]) -> str:
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    output = (result.stdout or result.stderr).strip().splitlines()
    return output[0] if output else f"exit={result.returncode}; no version text"


def build_run_manifest(inputs: list[Path], challenge_revisions: dict[str, str]) -> dict[str, Any]:
    version_commands = {
        "bcftools": ["bcftools", "--version"],
        "fastp": ["fastp", "--version"],
        "multiqc": ["multiqc", "--version"],
        "samtools": ["samtools", "--version"],
        "snakemake": ["snakemake", "--version"],
    }
    return {
        "schema_version": "1.0",
        "created_at": datetime.now(UTC).isoformat(),
        "challenge_revisions": challenge_revisions,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "tools": {
            name: tool_version(command) for name, command in sorted(version_commands.items())
        },
        "inputs": [
            {"name": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in inputs
        ],
    }


def write_manifest(manifest: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
