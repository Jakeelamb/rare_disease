#!/usr/bin/env python3
"""Exercise chromosome-scattered WhatsHap phasing and ordered BCFtools gather."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

BIN = Path(sys.executable).parent


def tool(name: str) -> str:
    path = BIN / name
    if not path.is_file():
        raise RuntimeError(f"required phasing executable is absent from Pixi environment: {path}")
    return str(path)


def run(command: list[str]) -> None:
    subprocess.run(command, check=True, capture_output=True, text=True)


def write_reference(path: Path) -> None:
    path.write_text(
        ">1\n" + "A" * 1000 + "\n>2\n" + "A" * 1000 + "\n>decoy\n" + "A" * 1000 + "\n",
        encoding="utf-8",
    )
    pysam.faidx(str(path))


def write_vcf(path: Path) -> Path:
    path.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##contig=<ID=2,length=1000>
##contig=<ID=decoy,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND01
1\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
1\t200\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
2\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
2\t200\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
decoy\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
""",
        encoding="utf-8",
    )
    compressed = path.with_suffix(".vcf.gz")
    pysam.tabix_compress(str(path), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def write_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "1", "LN": 1000},
            {"SN": "2", "LN": 1000},
            {"SN": "decoy", "LN": 1000},
        ],
        "RG": [{"ID": "synthetic", "SM": "PROBAND01"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for reference_id in (0, 1):
            for index, (left, right) in enumerate((("G", "A"), ("G", "A"), ("A", "G"), ("A", "G"))):
                sequence = list("A" * 250)
                sequence[49] = left
                sequence[149] = right
                read = pysam.AlignedSegment()
                read.query_name = f"synthetic_{reference_id}_{index}"
                read.query_sequence = "".join(sequence)
                read.flag = 0
                read.reference_id = reference_id
                read.reference_start = 50
                read.mapping_quality = 60
                read.cigar = ((0, 250),)
                read.query_qualities = pysam.qualitystring_to_array("I" * 250)
                read.set_tag("RG", "synthetic")
                bam.write(read)
    pysam.index(str(path))


def index_vcf(path: Path) -> None:
    run([tool("bcftools"), "index", "--tbi", "--output", f"{path}.tbi", str(path)])


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="mva-phasing-scatter-smoke-") as raw_directory:
        directory = Path(raw_directory)
        reference = directory / "reference.fa"
        input_vcf = write_vcf(directory / "input.vcf")
        bam = directory / "reads.bam"
        write_reference(reference)
        write_bam(bam)

        phased_shards: list[Path] = []
        for contig in ("1", "2"):
            split = directory / f"{contig}.input.vcf.gz"
            phased = directory / f"{contig}.phased.vcf.gz"
            run(
                [
                    tool("bcftools"),
                    "view",
                    "--regions",
                    contig,
                    "--output-type",
                    "z",
                    "--output",
                    str(split),
                    str(input_vcf),
                ]
            )
            index_vcf(split)
            run(
                [
                    tool("whatshap"),
                    "phase",
                    "--reference",
                    str(reference),
                    "--sample",
                    "PROBAND01",
                    "--output",
                    str(phased),
                    str(split),
                    str(bam),
                ]
            )
            index_vcf(phased)
            phased_shards.append(phased)

        nonprimary = directory / "nonprimary.unphased.vcf.gz"
        run(
            [
                tool("bcftools"),
                "view",
                "--targets",
                "^1,2",
                "--output-type",
                "z",
                "--output",
                str(nonprimary),
                str(input_vcf),
            ]
        )
        index_vcf(nonprimary)

        gathered = directory / "gathered.vcf.gz"
        run(
            [
                tool("bcftools"),
                "concat",
                "--output-type",
                "z",
                "--output",
                str(gathered),
                *(str(path) for path in (*phased_shards, nonprimary)),
            ]
        )
        index_vcf(gathered)

        with pysam.VariantFile(str(gathered)) as vcf:
            records = tuple(vcf)
        if [record.contig for record in records] != ["1", "1", "2", "2", "decoy"]:
            raise RuntimeError("gathered VCF did not preserve reference contig order")
        phased_primary = sum(
            record.samples["PROBAND01"].phased for record in records if record.contig in {"1", "2"}
        )
        if len(records) != 5 or phased_primary != 4:
            raise RuntimeError(
                f"unexpected gathered result: records={len(records)}, phased={phased_primary}"
            )
        if records[-1].samples["PROBAND01"].phased:
            raise RuntimeError("non-primary call was unexpectedly phased")

        print(
            json.dumps(
                {
                    "status": "PASS",
                    "records_preserved": len(records),
                    "primary_variants_phased": phased_primary,
                    "nonprimary_variants_preserved_unphased": 1,
                    "patient_data_used": False,
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
