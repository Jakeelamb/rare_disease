#!/usr/bin/env python3
"""Rebuild the reference named by the source VCF from primary public components."""

from __future__ import annotations

import argparse
import gzip
import json
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open(encoding="ascii")


def fasta_records(handle: TextIO) -> Iterator[tuple[str, list[str]]]:
    header: str | None = None
    sequence: list[str] = []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, sequence
            header = line.rstrip("\n")
            sequence = []
        else:
            if header is None:
                raise ValueError("sequence occurred before the first FASTA header")
            sequence.append(line.rstrip("\n"))
    if header is not None:
        yield header, sequence


def normalized_name(header: str) -> str:
    name = header[1:].split(maxsplit=1)[0]
    return name[3:] if name.startswith("chr") else name


def normalized_header(header: str) -> str:
    parts = header[1:].split(maxsplit=1)
    name = normalized_name(header)
    return f">{name} {parts[1]}" if len(parts) == 2 else f">{name}"


def sequence_length(sequence: list[str]) -> int:
    return sum(len(line.strip()) for line in sequence)


def build_reference(masked_primary: Path, combined: Path, output: Path) -> dict[str, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    primary_count = 0
    total_count = 0
    names: set[str] = set()

    with (
        open_text(masked_primary) as primary_handle,
        open_text(combined) as combined_handle,
        output.open("w", encoding="ascii") as destination,
    ):
        combined_records = fasta_records(combined_handle)
        for primary_header, primary_sequence in fasta_records(primary_handle):
            try:
                combined_header, combined_sequence = next(combined_records)
            except StopIteration as error:
                raise ValueError("combined reference ended before primary reference") from error
            primary_name = normalized_name(primary_header)
            combined_name = normalized_name(combined_header)
            if primary_name != combined_name:
                raise ValueError(f"primary order mismatch: {primary_name!r} != {combined_name!r}")
            if sequence_length(primary_sequence) != sequence_length(combined_sequence):
                raise ValueError(f"primary length mismatch for {primary_name}")
            if primary_name in names:
                raise ValueError(f"duplicate contig {primary_name}")
            names.add(primary_name)
            destination.write(normalized_header(primary_header) + "\n")
            destination.writelines(f"{line}\n" for line in primary_sequence)
            primary_count += 1
            total_count += 1

        for header, sequence in combined_records:
            name = normalized_name(header)
            if name in names:
                raise ValueError(f"duplicate contig {name}")
            names.add(name)
            destination.write(normalized_header(header) + "\n")
            destination.writelines(f"{line}\n" for line in sequence)
            total_count += 1

    return {"masked_primary_contigs": primary_count, "total_contigs": total_count}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masked-primary", type=Path, required=True)
    parser.add_argument("--combined", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = build_reference(args.masked_primary, args.combined, args.output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
