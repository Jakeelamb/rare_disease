#!/usr/bin/env python3
"""Exercise the production BAM command sequence on deterministic synthetic reads."""

from __future__ import annotations

import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = Path(sys.executable).parent


def tool(name: str) -> str:
    path = BIN / name
    if not path.is_file():
        raise RuntimeError(f"required WGS executable is absent from Pixi environment: {path}")
    return str(path)


def run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, **kwargs)  # type: ignore[arg-type]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for name, sequence in records:
        lines.extend((f"@{name}", sequence, "+", "I" * len(sequence)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def align_lane(reference: Path, r1: Path, r2: Path, output: Path, read_group: str) -> None:
    align = subprocess.Popen(
        [tool("bwa-mem2"), "mem", "-t", "2", "-R", read_group, str(reference), str(r1), str(r2)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    assert align.stdout is not None
    sort = subprocess.run(
        [tool("samtools"), "sort", "-n", "-@", "2", "-O", "BAM", "-o", str(output), "-"],
        stdin=align.stdout,
        check=True,
        capture_output=True,
    )
    align.stdout.close()
    stderr = align.stderr.read().decode() if align.stderr is not None else ""
    if align.wait() != 0:
        raise RuntimeError(f"bwa-mem2 failed: {stderr}")
    if sort.returncode:
        raise RuntimeError(sort.stderr.decode())


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="mva-wgs-smoke-") as raw_directory:
        directory = Path(raw_directory)
        random_generator = random.Random(20260825)
        reference_sequence = "".join(random_generator.choice("ACGT") for _ in range(20_000))
        reference = directory / "reference.fa"
        reference.write_text(">chrSynthetic\n" + reference_sequence + "\n", encoding="utf-8")
        run([tool("samtools"), "faidx", str(reference)], capture_output=True)
        run([tool("bwa-mem2"), "index", str(reference)], capture_output=True)

        lane_bams: list[Path] = []
        for lane_index, starts in enumerate(((1000, 2000), (3000, 4000)), start=1):
            r1_records: list[tuple[str, str]] = []
            r2_records: list[tuple[str, str]] = []
            for pair_index, start in enumerate(starts, start=1):
                name = f"synthetic_{lane_index}_{pair_index}"
                r1_records.append((f"{name}/1", reference_sequence[start : start + 150]))
                r2_records.append(
                    (f"{name}/2", reverse_complement(reference_sequence[start + 350 : start + 500]))
                )
            r1 = directory / f"lane{lane_index}_R1.fastq"
            r2 = directory / f"lane{lane_index}_R2.fastq"
            write_fastq(r1, r1_records)
            write_fastq(r2, r2_records)
            lane_bam = directory / f"lane{lane_index}.name.bam"
            align_lane(
                reference,
                r1,
                r2,
                lane_bam,
                f"@RG\tID:L{lane_index}\tSM:SYNTHETIC\tLB:WGS\tPL:ILLUMINA",
            )
            lane_bams.append(lane_bam)

        merged = directory / "merged.name.bam"
        run(
            [
                tool("samtools"),
                "merge",
                "-n",
                "-@",
                "2",
                "-O",
                "BAM",
                "-o",
                str(merged),
                *(str(path) for path in lane_bams),
            ],
            capture_output=True,
        )
        fixmate = subprocess.Popen(
            [tool("samtools"), "fixmate", "-m", "-@", "2", str(merged), "-"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert fixmate.stdout is not None
        coordinate = directory / "coordinate.bam"
        subprocess.run(
            [
                tool("samtools"),
                "sort",
                "-@",
                "2",
                "-O",
                "BAM",
                "-o",
                str(coordinate),
                "-",
            ],
            stdin=fixmate.stdout,
            check=True,
            capture_output=True,
        )
        fixmate.stdout.close()
        fixmate_stderr = fixmate.stderr.read().decode() if fixmate.stderr is not None else ""
        if fixmate.wait() != 0:
            raise RuntimeError(f"samtools fixmate failed: {fixmate_stderr}")

        marked = directory / "marked.bam"
        run(
            [tool("samtools"), "markdup", "-s", "-@", "2", str(coordinate), str(marked)],
            capture_output=True,
        )
        run([tool("samtools"), "index", "-@", "2", str(marked)], capture_output=True)
        run([tool("samtools"), "quickcheck", "-v", str(marked)], capture_output=True)
        flagstat = run(
            [tool("samtools"), "flagstat", "-@", "2", str(marked)], capture_output=True
        ).stdout
        if "8 + 0 in total" not in flagstat:
            raise RuntimeError(f"unexpected flagstat output:\n{flagstat}")
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "synthetic_reads": 8,
                    "samtools": run(
                        [tool("samtools"), "--version"], capture_output=True
                    ).stdout.splitlines()[0],
                    "bwa_mem2": "2.3 package; executable capability tested",
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
