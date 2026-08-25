"""Privacy-minimizing aggregate inspection of the distributed VCF."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from cyvcf2 import VCF  # type: ignore[import-untyped]


def inspect_vcf(path: Path) -> dict[str, Any]:
    vcf = VCF(str(path))
    samples = tuple(vcf.samples)
    contigs = tuple(vcf.seqnames)
    filters: Counter[str] = Counter()
    genotypes: Counter[str] = Counter()
    records = 0
    snvs = 0
    indels = 0
    multiallelic = 0

    for record in vcf:
        records += 1
        filters[record.FILTER or "PASS"] += 1
        multiallelic += int(len(record.ALT) > 1)
        if len(record.REF) == 1 and all(len(alt) == 1 for alt in record.ALT):
            snvs += 1
        else:
            indels += 1
        if record.genotypes:
            left, right, phased = record.genotypes[0]
            separator = "|" if phased else "/"
            genotypes[f"{left}{separator}{right}"] += 1
    vcf.close()
    return {
        "path": str(path),
        "sample_count": len(samples),
        "samples": samples,
        "contig_count": len(contigs),
        "primary_contigs": tuple(
            contig
            for contig in contigs
            if contig in {str(i) for i in range(1, 23)} | {"X", "Y", "M"}
        ),
        "uses_chr_prefix": any(contig.startswith("chr") for contig in contigs),
        "record_count": records,
        "snv_count": snvs,
        "indel_or_complex_count": indels,
        "multiallelic_count": multiallelic,
        "filter_counts": dict(filters.most_common()),
        "genotype_counts": dict(genotypes.most_common()),
    }
