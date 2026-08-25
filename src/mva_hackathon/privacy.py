"""Fail-closed checks for files intended to enter the public repository."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PrivacyFinding:
    rule: str
    path: str
    detail: str


@dataclass(frozen=True)
class PrivacyAuditResult:
    findings: tuple[PrivacyFinding, ...]

    @property
    def ok(self) -> bool:
        return not self.findings


RESTRICTED_PREFIXES = ("data/", "results/private/", "results/runs/", "work/")
RESTRICTED_SUFFIXES = (
    ".bam",
    ".bcf",
    ".cram",
    ".docx",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".g.vcf",
    ".g.vcf.gz",
    ".vcf",
    ".vcf.gz",
)
SECRET_PATTERNS = {
    "secret.huggingface_token": re.compile(rb"hf_[A-Za-z0-9]{20,}"),
    "secret.github_token": re.compile(rb"gh[pousr]_[A-Za-z0-9]{30,}"),
    "secret.aws_access_key": re.compile(rb"AKIA[0-9A-Z]{16}"),
}
MAX_PUBLIC_FILE_BYTES = 5 * 1024 * 1024


def public_candidate_files(root: Path) -> list[Path]:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        detail = result.stderr.decode().strip()
        raise ValueError(f"cannot enumerate public Git candidates: {detail}")
    return [root / item.decode() for item in result.stdout.split(b"\0") if item]


def audit_public_tree(
    root: Path,
    *,
    files: list[Path] | None = None,
    max_scan_bytes: int = 10 * 1024 * 1024,
    max_public_file_bytes: int = MAX_PUBLIC_FILE_BYTES,
) -> PrivacyAuditResult:
    root = root.resolve()
    candidates = files if files is not None else public_candidate_files(root)
    findings: list[PrivacyFinding] = []
    for path in candidates:
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(root).as_posix()
        except ValueError:
            relative = resolved.name
        lowered = relative.lower()
        if lowered.startswith(RESTRICTED_PREFIXES):
            findings.append(
                PrivacyFinding(
                    rule="restricted.path",
                    path=relative,
                    detail="restricted patient-derived path is visible to Git",
                )
            )
        if lowered.endswith(RESTRICTED_SUFFIXES):
            findings.append(
                PrivacyFinding(
                    rule="restricted.genomic_extension",
                    path=relative,
                    detail="patient-capable binary/genomic artifact must remain ignored",
                )
            )
        if not resolved.is_file():
            continue
        size = resolved.stat().st_size
        if size > max_public_file_bytes:
            findings.append(
                PrivacyFinding(
                    rule="restricted.oversized_file",
                    path=relative,
                    detail=(
                        f"public Git candidate is {size} bytes; "
                        f"limit is {max_public_file_bytes} bytes"
                    ),
                )
            )
        if size > max_scan_bytes:
            continue
        content = resolved.read_bytes()
        for rule, pattern in SECRET_PATTERNS.items():
            if pattern.search(content):
                findings.append(
                    PrivacyFinding(
                        rule=rule,
                        path=relative,
                        detail="high-confidence credential pattern detected",
                    )
                )
    return PrivacyAuditResult(
        findings=tuple(sorted(findings, key=lambda item: (item.path, item.rule)))
    )
