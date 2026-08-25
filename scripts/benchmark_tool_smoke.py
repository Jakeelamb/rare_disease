#!/usr/bin/env python3
"""Exercise the production vcfeval adapter on a tiny public synthetic truth set."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

from run_vcfeval import run_vcfeval


def write_vcf(path: Path, rows: list[str]) -> None:
    header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
    subprocess.run(["bgzip", "--force", str(path)], check=True)
    subprocess.run(["tabix", "--preset", "vcf", f"{path}.gz"], check=True)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="mva-benchmark-smoke-") as raw_directory:
        directory = Path(raw_directory)
        reference = directory / "reference.fa"
        reference.write_text(">chr1\n" + "ACGT" * 250 + "\n", encoding="utf-8")
        truth = directory / "truth.vcf"
        calls = directory / "calls.vcf"
        write_vcf(truth, ["chr1\t101\t.\tA\tG\t60\tPASS\t.\tGT:GQ\t0/1:60"])
        write_vcf(
            calls,
            [
                "chr1\t101\t.\tA\tG\t60\tPASS\t.\tGT:GQ\t0/1:60",
                "chr1\t201\t.\tA\tT\t60\tPASS\t.\tGT:GQ\t0/1:60",
            ],
        )
        metrics_path = run_vcfeval(
            truth=Path(f"{truth}.gz"),
            calls=Path(f"{calls}.gz"),
            reference=reference,
            output=directory / "result",
            truth_name="synthetic-GRCh38-like",
            callset_name="synthetic-caller",
        )
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        observed = (
            metrics["true_positive"],
            metrics["false_positive"],
            metrics["false_negative"],
        )
        if observed != (1, 1, 0):
            raise RuntimeError(f"unexpected vcfeval smoke counts: {observed}")
        print(json.dumps({"status": "PASS", "tp": 1, "fp": 1, "fn": 0}, indent=2))


if __name__ == "__main__":
    main()
