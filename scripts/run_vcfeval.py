#!/usr/bin/env python3
"""Run one public truth-set vcfeval comparison and emit standardized metrics."""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
import time
from pathlib import Path


def command_output(command: list[str]) -> str:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return (result.stdout or result.stderr).strip().splitlines()[0]


def record_count(path: Path) -> int:
    result = subprocess.run(
        ["bcftools", "view", "--no-header", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return sum(bool(line) for line in result.stdout.splitlines())


def run_vcfeval(
    *,
    truth: Path,
    calls: Path,
    reference: Path,
    output: Path,
    truth_name: str,
    callset_name: str,
    evaluation_regions: Path | None = None,
) -> Path:
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"vcfeval output directory must be absent or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    sdf = output / "reference.sdf"
    subprocess.run(["rtg", "format", "--output", str(sdf), str(reference)], check=True)
    command = [
        "rtg",
        "vcfeval",
        "--baseline",
        str(truth),
        "--calls",
        str(calls),
        "--template",
        str(sdf),
        "--output",
        str(output / "vcfeval"),
        "--squash-ploidy",
    ]
    if evaluation_regions is not None:
        command.extend(["--evaluation-regions", str(evaluation_regions)])
    started = time.perf_counter()
    subprocess.run(command, check=True)
    runtime = time.perf_counter() - started
    peak_rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    result_dir = output / "vcfeval"
    metrics = {
        "callset": callset_name,
        "truth_set": truth_name,
        "engine": command_output(["rtg", "version"]),
        "true_positive": record_count(result_dir / "tp-baseline.vcf.gz"),
        "false_positive": record_count(result_dir / "fp.vcf.gz"),
        "false_negative": record_count(result_dir / "fn.vcf.gz"),
        "runtime_seconds": runtime,
        "peak_rss_mb": max(peak_rss_mb, 0.001),
        "command": command,
    }
    metrics_path = output / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    return metrics_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--calls", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--truth-name", required=True)
    parser.add_argument("--callset-name", required=True)
    parser.add_argument("--evaluation-regions", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    path = run_vcfeval(
        truth=args.truth,
        calls=args.calls,
        reference=args.reference,
        output=args.output,
        truth_name=args.truth_name,
        callset_name=args.callset_name,
        evaluation_regions=args.evaluation_regions,
    )
    print(path.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
