"""Candidate-level BAM evidence without exposing read names or sequences."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pysam
from pydantic import BaseModel, Field

from .models import RankedCase, VariantKey


class AlleleReadEvidence(BaseModel):
    variant: VariantKey
    usable_depth: int = Field(ge=0)
    ref_reads: int = Field(ge=0)
    alt_reads: int = Field(ge=0)
    other_reads: int = Field(ge=0)
    filtered_reads: int = Field(ge=0)
    alt_forward_reads: int = Field(ge=0)
    alt_reverse_reads: int = Field(ge=0)
    alt_proper_pair_reads: int = Field(default=0, ge=0)
    alt_soft_clipped_reads: int = Field(default=0, ge=0)
    mean_alt_mapping_quality: float | None = None
    mean_alt_base_quality: float | None = None
    mean_alt_distance_from_read_end: float | None = None

    @property
    def allele_balance(self) -> float | None:
        informative = self.ref_reads + self.alt_reads
        return self.alt_reads / informative if informative else None


class PairReadEvidence(BaseModel):
    rank: int | None = Field(default=None, ge=1)
    gene: str | None = None
    left: VariantKey
    right: VariantKey
    left_alt_observed: bool = True
    right_alt_observed: bool = True
    shared_fragments: int = Field(ge=0)
    alt_alt_fragments: int = Field(ge=0)
    alt_ref_fragments: int = Field(ge=0)
    ref_alt_fragments: int = Field(ge=0)
    ref_ref_fragments: int = Field(ge=0)
    other_fragments: int = Field(ge=0)
    min_phase_support_fragments: int = Field(default=2, ge=1)

    @property
    def interpretation(self) -> str:
        if not self.left_alt_observed or not self.right_alt_observed:
            return "insufficient_exact_alt_support"
        trans = self.alt_ref_fragments + self.ref_alt_fragments
        if self.alt_alt_fragments and trans:
            return "conflicting_read_phase"
        if self.alt_alt_fragments >= self.min_phase_support_fragments:
            return "read_backed_cis"
        if trans >= self.min_phase_support_fragments:
            return "read_backed_trans"
        if self.alt_alt_fragments or trans:
            return "insufficient_phase_support"
        if self.shared_fragments:
            return "no_informative_spanning_fragments"
        return "no_spanning_fragments"


class CandidateReadAudit(BaseModel):
    bam: str
    min_mapping_quality: int
    min_base_quality: int
    min_phase_support_fragments: int = Field(default=2, ge=1)
    alleles: tuple[AlleleReadEvidence, ...]
    pairs: tuple[PairReadEvidence, ...]


def _alignment_contig(bam: pysam.AlignmentFile, key: VariantKey) -> str:
    candidates = (key.chrom, key.challenge_chrom)
    for candidate in candidates:
        if candidate in bam.references:
            return candidate
    raise ValueError(f"BAM has no contig for {key.challenge_chrom}")


def _observed_allele(
    pileup_read: pysam.PileupRead, key: VariantKey, *, min_base_quality: int
) -> str:
    query_position = pileup_read.query_position
    alignment = pileup_read.alignment
    sequence = alignment.query_sequence
    qualities = alignment.query_qualities
    if pileup_read.is_refskip or query_position is None or sequence is None:
        return "other"
    if qualities is None or qualities[query_position] < min_base_quality:
        return "filtered"

    ref = key.ref
    alt = key.alt
    if len(ref) == len(alt):
        observed = sequence[query_position : query_position + len(ref)].upper()
        if observed == ref:
            return "ref"
        if observed == alt:
            return "alt"
        return "other"

    if len(alt) > len(ref) and alt.startswith(ref):
        inserted = alt[len(ref) :]
        observed = sequence[query_position + 1 : query_position + 1 + len(inserted)].upper()
        if pileup_read.indel == len(inserted) and observed == inserted:
            return "alt"
        return "ref" if pileup_read.indel == 0 else "other"

    if len(ref) > len(alt) and ref.startswith(alt):
        deleted = len(ref) - len(alt)
        if pileup_read.indel == -deleted:
            return "alt"
        return "ref" if pileup_read.indel == 0 else "other"

    return "other"


def _inspect_allele(
    bam: pysam.AlignmentFile,
    key: VariantKey,
    *,
    min_mapping_quality: int,
    min_base_quality: int,
) -> tuple[AlleleReadEvidence, dict[str, str]]:
    counts = {label: 0 for label in ("ref", "alt", "other", "filtered")}
    alt_mapping_qualities: list[int] = []
    alt_base_qualities: list[int] = []
    alt_distances_from_read_end: list[int] = []
    alt_forward = 0
    alt_reverse = 0
    alt_proper_pair = 0
    alt_soft_clipped = 0
    fragment_calls: dict[str, str] = {}
    contig = _alignment_contig(bam, key)
    for column in bam.pileup(
        contig,
        key.pos - 1,
        key.pos,
        truncate=True,
        stepper="all",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=100_000,
    ):
        if column.reference_pos != key.pos - 1:
            continue
        for pileup_read in column.pileups:
            alignment = pileup_read.alignment
            if (
                alignment.is_unmapped
                or alignment.is_secondary
                or alignment.is_supplementary
                or alignment.is_duplicate
                or alignment.is_qcfail
                or alignment.mapping_quality < min_mapping_quality
            ):
                counts["filtered"] += 1
                continue
            call = _observed_allele(pileup_read, key, min_base_quality=min_base_quality)
            counts[call] += 1
            query_name = alignment.query_name
            if call in {"ref", "alt"} and query_name is not None:
                previous = fragment_calls.get(query_name)
                fragment_calls[query_name] = call if previous in {None, call} else "other"
            if call == "alt":
                alt_reverse += int(alignment.is_reverse)
                alt_forward += int(not alignment.is_reverse)
                alt_proper_pair += int(alignment.is_proper_pair)
                alt_soft_clipped += int(
                    any(operation == 4 for operation, _ in (alignment.cigartuples or ()))
                )
                alt_mapping_qualities.append(alignment.mapping_quality)
                query_position = pileup_read.query_position
                if alignment.query_qualities is not None and query_position is not None:
                    alt_base_qualities.append(alignment.query_qualities[query_position])
                if alignment.query_length is not None and query_position is not None:
                    alt_distances_from_read_end.append(
                        min(query_position, alignment.query_length - query_position - 1)
                    )

    evidence = AlleleReadEvidence(
        variant=key,
        usable_depth=counts["ref"] + counts["alt"] + counts["other"],
        ref_reads=counts["ref"],
        alt_reads=counts["alt"],
        other_reads=counts["other"],
        filtered_reads=counts["filtered"],
        alt_forward_reads=alt_forward,
        alt_reverse_reads=alt_reverse,
        alt_proper_pair_reads=alt_proper_pair,
        alt_soft_clipped_reads=alt_soft_clipped,
        mean_alt_mapping_quality=(
            round(sum(alt_mapping_qualities) / len(alt_mapping_qualities), 3)
            if alt_mapping_qualities
            else None
        ),
        mean_alt_base_quality=(
            round(sum(alt_base_qualities) / len(alt_base_qualities), 3)
            if alt_base_qualities
            else None
        ),
        mean_alt_distance_from_read_end=(
            round(sum(alt_distances_from_read_end) / len(alt_distances_from_read_end), 3)
            if alt_distances_from_read_end
            else None
        ),
    )
    return evidence, fragment_calls


def _inspect_pair(
    left: VariantKey,
    right: VariantKey,
    calls: Mapping[VariantKey, Mapping[str, str]],
    *,
    rank: int | None = None,
    gene: str | None = None,
    min_phase_support_fragments: int = 2,
) -> PairReadEvidence:
    shared = set(calls[left]) & set(calls[right])
    combinations = {label: 0 for label in ("alt_alt", "alt_ref", "ref_alt", "ref_ref", "other")}
    for fragment in shared:
        combination = f"{calls[left][fragment]}_{calls[right][fragment]}"
        combinations[combination if combination in combinations else "other"] += 1
    return PairReadEvidence(
        rank=rank,
        gene=gene,
        left=left,
        right=right,
        left_alt_observed=any(call == "alt" for call in calls[left].values()),
        right_alt_observed=any(call == "alt" for call in calls[right].values()),
        shared_fragments=len(shared),
        alt_alt_fragments=combinations["alt_alt"],
        alt_ref_fragments=combinations["alt_ref"],
        ref_alt_fragments=combinations["ref_alt"],
        ref_ref_fragments=combinations["ref_ref"],
        other_fragments=combinations["other"],
        min_phase_support_fragments=min_phase_support_fragments,
    )


def inspect_ranked_bam(
    ranked: RankedCase,
    bam_path: Path,
    *,
    limit: int = 10,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_phase_support_fragments: int = 2,
) -> CandidateReadAudit:
    candidates = ranked.candidates[:limit]
    keys = tuple(dict.fromkeys(variant.key for item in candidates for variant in item.variants))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        inspected = [
            _inspect_allele(
                bam,
                key,
                min_mapping_quality=min_mapping_quality,
                min_base_quality=min_base_quality,
            )
            for key in keys
        ]
    evidence = tuple(item[0] for item in inspected)
    calls = {key: inspected[index][1] for index, key in enumerate(keys)}
    pairs = tuple(
        _inspect_pair(
            item.variants[0].key,
            item.variants[1].key,
            calls,
            rank=rank,
            gene=item.gene,
            min_phase_support_fragments=min_phase_support_fragments,
        )
        for rank, item in enumerate(candidates, start=1)
        if len(item.variants) == 2
    )
    return CandidateReadAudit(
        bam=str(bam_path),
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_phase_support_fragments=min_phase_support_fragments,
        alleles=evidence,
        pairs=pairs,
    )


def render_candidate_read_audit(audit: CandidateReadAudit) -> str:
    lines = [
        "# Private candidate read-evidence review",
        "",
        f"- Minimum mapping quality: {audit.min_mapping_quality}",
        f"- Minimum base quality: {audit.min_base_quality}",
        "- Minimum concordant informative fragments for a phase call: "
        f"{audit.min_phase_support_fragments}",
        "",
        "| Allele | Depth | Ref | Alt | AB | Alt F/R | Proper pair | Soft clipped | "
        "Mean alt MQ | Mean alt BQ | Mean end distance |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for allele in audit.alleles:
        balance = f"{allele.allele_balance:.3f}" if allele.allele_balance is not None else "n/a"
        mean_mapping_quality = (
            allele.mean_alt_mapping_quality
            if allele.mean_alt_mapping_quality is not None
            else "n/a"
        )
        mean_base_quality = (
            allele.mean_alt_base_quality if allele.mean_alt_base_quality is not None else "n/a"
        )
        mean_end_distance = (
            allele.mean_alt_distance_from_read_end
            if allele.mean_alt_distance_from_read_end is not None
            else "n/a"
        )
        lines.append(
            f"| `{allele.variant.label}` | {allele.usable_depth} | {allele.ref_reads} | "
            f"{allele.alt_reads} | {balance} | "
            f"{allele.alt_forward_reads}/{allele.alt_reverse_reads} | "
            f"{allele.alt_proper_pair_reads} | {allele.alt_soft_clipped_reads} | "
            f"{mean_mapping_quality} | {mean_base_quality} | {mean_end_distance} |"
        )
    lines.extend(
        [
            "",
            "| Rank | Gene | Pair | Exact alt L/R | Shared fragments | Alt/alt | Alt/ref | "
            "Ref/alt | Ref/ref | Other | Interpretation |",
            "|---:|---|---|---|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for pair in audit.pairs:
        lines.append(
            f"| {pair.rank or 'n/a'} | {pair.gene or 'n/a'} | "
            f"`{pair.left.label}` + `{pair.right.label}` | "
            f"{int(pair.left_alt_observed)}/{int(pair.right_alt_observed)} | "
            f"{pair.shared_fragments} | "
            f"{pair.alt_alt_fragments} | {pair.alt_ref_fragments} | "
            f"{pair.ref_alt_fragments} | {pair.ref_ref_fragments} | "
            f"{pair.other_fragments} | `{pair.interpretation}` |"
        )
    return "\n".join(lines) + "\n"


def write_candidate_read_audit(
    audit: CandidateReadAudit, *, json_output: Path, review_output: Path
) -> None:
    json_output.parent.mkdir(parents=True, exist_ok=True)
    review_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(audit.model_dump_json(indent=2) + "\n", encoding="utf-8")
    review_output.write_text(render_candidate_read_audit(audit), encoding="utf-8")
