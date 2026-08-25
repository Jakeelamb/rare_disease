#!/usr/bin/env python3
"""Run the lightweight workflow stages on synthetic local data."""

from __future__ import annotations

import gzip
import json
import subprocess
import tempfile
from pathlib import Path

import yaml
from docx import Document


def write_phenotype(path: Path) -> None:
    document = Document()
    table = document.add_table(rows=2, cols=4)
    headers = ("Clinical Feature", "HPO Term", "HPO ID", "Presentation / Notes")
    for cell, value in zip(table.rows[0].cells, headers, strict=True):
        cell.text = value
    values = ("Synthetic seizure", "Seizure", "HP:0001250", "public synthetic fixture")
    for cell, value in zip(table.rows[1].cells, values, strict=True):
        cell.text = value
    document.save(path)


def write_fastq(path: Path, mate: int) -> None:
    sequence = "ACGT" * 25 if mate == 1 else "TGCA" * 25
    content = f"@synthetic/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(content)


def run_snakemake(workspace: Path, config_path: Path, target: str) -> None:
    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "4",
            target,
        ],
        cwd=workspace,
        check=True,
    )


def main() -> None:
    workspace = Path.cwd()
    artifacts = workspace / ".artifacts"
    artifacts.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="workflow-smoke-", dir=artifacts) as raw:
        root = Path(raw)
        data = root / "data"
        data.mkdir()
        phenotype = data / "phenotype.docx"
        write_phenotype(phenotype)
        r1 = data / "synthetic_L001_R1_001.fastq.gz"
        r2 = data / "synthetic_L001_R2_001.fastq.gz"
        write_fastq(r1, 1)
        write_fastq(r2, 2)
        vcf = data / "proband.vcf"
        vcf.write_text(
            """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
1\t101\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
""",
            encoding="utf-8",
        )
        subprocess.run(["bgzip", "--force", str(vcf)], check=True)
        compressed_vcf = Path(f"{vcf}.gz")
        subprocess.run(["tabix", "--preset", "vcf", str(compressed_vcf)], check=True)

        config = yaml.safe_load(Path("config/workflow.yaml").read_text(encoding="utf-8"))
        relative_root = root.relative_to(workspace).as_posix()
        config["inputs"] = {
            "phenotype_glob": f"{relative_root}/data/*.docx",
            "phenotype_placeholder": f"{relative_root}/data/phenotype.docx",
            "vcf_glob": f"{relative_root}/data/*.vcf.gz",
            "vcf_placeholder": f"{relative_root}/data/proband.vcf.gz",
            "fastq_glob": f"{relative_root}/data/*.fastq.gz",
        }
        config["outputs"]["private_root"] = f"{relative_root}/work"
        config_path = root / "workflow.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

        run_snakemake(workspace, config_path, "intake")
        run_snakemake(workspace, config_path, "phenotype_extraction")

        work = root / "work"
        expected = [
            work / "qc" / "multiqc" / "multiqc_report.html",
            work / "qc" / "vcf_summary.json",
            work / "phenotype" / "phenotype.json",
            work / "phenotype" / "phenotype.review.md",
        ]
        missing = [str(path) for path in expected if not path.is_file()]
        if missing:
            raise RuntimeError(f"workflow smoke outputs are missing: {missing}")
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "stages": ["intake", "phenotype_extraction"],
                    "patient_data_used": False,
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
