"""Exact, local evidence integration for manual leading-candidate adjudication."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Any, Literal

import pysam
from pydantic import BaseModel, Field

from .models import CandidateHypothesis, RankedCase, VariantKey
from .read_evidence import AlleleReadEvidence, CandidateReadAudit, PairReadEvidence

type RecallStatus = Literal[
    "both_exact_alleles_recalled",
    "one_exact_allele_recalled",
    "no_exact_alleles_recalled",
    "single_exact_allele_recalled",
    "single_exact_allele_not_recalled",
]
type RecallPairPhase = Literal[
    "read_backed_trans",
    "read_backed_cis",
    "no_shared_phase_block",
    "not_applicable",
]


class ClinVarMatch(BaseModel):
    """One exact GRCh38 allele match from the local ClinVar variant summary."""

    variant: VariantKey
    allele_id: str
    variation_id: str
    name: str
    gene: str
    clinical_significance: str
    review_status: str
    number_submitters: int | None = Field(default=None, ge=0)
    last_evaluated: str | None = None
    rs_id: str | None = None
    rcv_accessions: tuple[str, ...] = ()
    phenotypes: tuple[str, ...] = ()


class ClinVarCandidateQuery(BaseModel):
    """Versioned exact-match results, including explicit no-match variants."""

    assembly: Literal["GRCh38"] = "GRCh38"
    release: str
    archive: str
    queried_variants: tuple[VariantKey, ...]
    matches: tuple[ClinVarMatch, ...]


class VcfCallObservation(BaseModel):
    """An exact allele observation from one indexed single-sample VCF."""

    caller: str
    variant: VariantKey
    present: bool
    sample: str | None = None
    gt: str | None = None
    dp: int | None = Field(default=None, ge=0)
    ad_ref: int | None = Field(default=None, ge=0)
    ad_alt: int | None = Field(default=None, ge=0)
    gq: float | None = Field(default=None, ge=0)
    qual: float | None = None
    filter: str | None = None
    phase_set: str | None = None
    phased: bool = False
    alt_haplotype: int | None = Field(default=None, ge=0, le=1)

    @property
    def allele_balance(self) -> float | None:
        if self.ad_ref is None or self.ad_alt is None:
            return None
        total = self.ad_ref + self.ad_alt
        return self.ad_alt / total if total else None


class LeadingCandidateAdjudication(BaseModel):
    """All evidence needed for a human decision on one exact hypothesis."""

    leading_candidate: CandidateHypothesis
    lane_ranks: dict[str, int | None]
    recall_calls: tuple[VcfCallObservation, ...]
    recall_status: RecallStatus
    recall_pair_phase: RecallPairPhase
    allele_read_evidence: tuple[AlleleReadEvidence, ...]
    pair_read_evidence: PairReadEvidence | None = None
    clinvar_release: str
    clinvar_matches: tuple[ClinVarMatch, ...]


def _candidate_identity(candidate: CandidateHypothesis) -> tuple[str, ...]:
    return tuple(sorted(variant.key.label for variant in candidate.variants))


def leading_variant_keys(
    ranked: RankedCase, *, hypothesis_limit: int = 3
) -> tuple[VariantKey, ...]:
    """Return exact variants from the leading hypotheses in stable first-seen order."""

    return tuple(
        dict.fromkeys(
            variant.key
            for candidate in ranked.candidates[:hypothesis_limit]
            for variant in candidate.variants
        )
    )


def _optional_int(value: str | None) -> int | None:
    if value is None or value in {"", "-", "na"}:
        return None
    return int(value)


def query_clinvar_archive(
    variants: tuple[VariantKey, ...], archive: Path, *, release: str
) -> ClinVarCandidateQuery:
    """Stream the public archive once and retain only exact GRCh38 allele matches."""

    targets = {
        (variant.chrom, variant.pos, variant.ref, variant.alt): variant for variant in variants
    }
    matches: list[ClinVarMatch] = []
    with gzip.open(archive, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["Assembly"] != "GRCh38":
                continue
            position = _optional_int(row.get("PositionVCF"))
            ref = row.get("ReferenceAlleleVCF")
            alt = row.get("AlternateAlleleVCF")
            chromosome = row.get("Chromosome")
            if position is None or not chromosome or not ref or not alt:
                continue
            try:
                key = VariantKey(chrom=chromosome, pos=position, ref=ref, alt=alt)
            except ValueError:
                continue
            target = targets.get((key.chrom, key.pos, key.ref, key.alt))
            if target is None:
                continue
            submitters = _optional_int(row.get("NumberSubmitters"))
            rs_number = row.get("RS# (dbSNP)")
            matches.append(
                ClinVarMatch(
                    variant=target,
                    allele_id=row["#AlleleID"],
                    variation_id=row["VariationID"],
                    name=row["Name"],
                    gene=row["GeneSymbol"],
                    clinical_significance=row["ClinicalSignificance"],
                    review_status=row["ReviewStatus"],
                    number_submitters=submitters,
                    last_evaluated=row.get("LastEvaluated") or None,
                    rs_id=(f"rs{rs_number}" if rs_number not in {None, "", "-"} else None),
                    rcv_accessions=tuple(
                        item for item in row.get("RCVaccession", "").split("|") if item != "-"
                    ),
                    phenotypes=tuple(
                        item for item in row.get("PhenotypeList", "").split("|") if item != "-"
                    ),
                )
            )
    matches.sort(key=lambda match: (match.variant.label, match.variation_id))
    return ClinVarCandidateQuery(
        release=release,
        archive=str(archive),
        queried_variants=variants,
        matches=tuple(matches),
    )


def write_clinvar_query(query: ClinVarCandidateQuery, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(query.model_dump_json(indent=2) + "\n", encoding="utf-8")


def _record_value(sample: Any, name: str) -> Any:
    return sample.get(name)


def _nonnegative_int(value: Any) -> int | None:
    if value is None:
        return None
    number = int(value)
    return number if number >= 0 else None


def _nonnegative_float(value: Any) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if number >= 0 else None


def _absent_call(caller: str, key: VariantKey) -> VcfCallObservation:
    return VcfCallObservation(caller=caller, variant=key, present=False)


def _call_from_record(
    caller: str,
    key: VariantKey,
    record: Any,
    sample_name: str,
) -> VcfCallObservation:
    alts = tuple(record.alts or ())
    alt_index = alts.index(key.alt) + 1
    sample = record.samples[sample_name]
    raw_gt = _record_value(sample, "GT")
    gt_values = tuple(raw_gt) if isinstance(raw_gt, tuple) else ()
    separator = "|" if sample.phased else "/"
    gt = separator.join("." if value is None else str(value) for value in gt_values) or None
    ad = _record_value(sample, "AD")
    ad_values = tuple(ad) if isinstance(ad, tuple) else ()
    phase_set = _record_value(sample, "PS") or _record_value(sample, "PID")
    alt_haplotype: int | None = None
    if sample.phased and len(gt_values) == 2 and gt_values.count(alt_index) == 1:
        alt_haplotype = gt_values.index(alt_index)
    filters = tuple(record.filter.keys())
    return VcfCallObservation(
        caller=caller,
        variant=key,
        present=alt_index in gt_values,
        sample=sample_name,
        gt=gt,
        dp=_nonnegative_int(_record_value(sample, "DP")),
        ad_ref=_nonnegative_int(ad_values[0]) if ad_values else None,
        ad_alt=_nonnegative_int(ad_values[alt_index]) if len(ad_values) > alt_index else None,
        gq=_nonnegative_float(_record_value(sample, "GQ")),
        qual=float(record.qual) if record.qual is not None else None,
        filter=";".join(filters) if filters else ".",
        phase_set=str(phase_set) if phase_set not in {None, ""} else None,
        phased=sample.phased,
        alt_haplotype=alt_haplotype,
    )


def inspect_exact_vcf_calls(
    vcf_path: Path,
    variants: tuple[VariantKey, ...],
    *,
    caller: str,
    sample_name: str | None = None,
) -> tuple[VcfCallObservation, ...]:
    """Inspect exact alleles without treating nearby or equivalent-looking rows as matches."""

    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = tuple(vcf.header.samples)
        if sample_name is None:
            if len(samples) != 1:
                raise ValueError("VCF must have exactly one sample unless sample_name is supplied")
            sample_name = samples[0]
        elif sample_name not in samples:
            raise ValueError(f"sample {sample_name} is absent from VCF")

        observations: list[VcfCallObservation] = []
        contigs = set(vcf.header.contigs)
        for key in variants:
            contig = next(
                (
                    candidate
                    for candidate in (key.chrom, key.challenge_chrom)
                    if candidate in contigs
                ),
                None,
            )
            if contig is None:
                observations.append(_absent_call(caller, key))
                continue
            exact_records = [
                record
                for record in vcf.fetch(contig, key.pos - 1, key.pos)
                if record.pos == key.pos
                and record.ref == key.ref
                and key.alt in (record.alts or ())
            ]
            if len(exact_records) > 1:
                raise ValueError(f"multiple exact VCF records found for {key.label}")
            observations.append(
                _call_from_record(caller, key, exact_records[0], sample_name)
                if exact_records
                else _absent_call(caller, key)
            )
    return tuple(observations)


def _recall_status(calls: tuple[VcfCallObservation, ...]) -> RecallStatus:
    present = sum(call.present for call in calls)
    if len(calls) == 1:
        return "single_exact_allele_recalled" if present else "single_exact_allele_not_recalled"
    if present == len(calls):
        return "both_exact_alleles_recalled"
    if present:
        return "one_exact_allele_recalled"
    return "no_exact_alleles_recalled"


def _recall_pair_phase(calls: tuple[VcfCallObservation, ...]) -> RecallPairPhase:
    if len(calls) != 2:
        return "not_applicable"
    left, right = calls
    if (
        not left.present
        or not right.present
        or left.phase_set is None
        or right.phase_set is None
        or left.phase_set != right.phase_set
        or left.alt_haplotype is None
        or right.alt_haplotype is None
    ):
        return "no_shared_phase_block"
    return "read_backed_cis" if left.alt_haplotype == right.alt_haplotype else "read_backed_trans"


def build_leading_adjudication(
    focused: RankedCase,
    lanes: dict[str, RankedCase],
    recall_calls: tuple[VcfCallObservation, ...],
    read_audit: CandidateReadAudit,
    clinvar: ClinVarCandidateQuery,
) -> LeadingCandidateAdjudication:
    if not focused.candidates:
        raise ValueError("focused ranking has no candidate to adjudicate")
    leading = focused.candidates[0]
    identity = _candidate_identity(leading)
    expected_keys = tuple(variant.key for variant in leading.variants)
    if tuple(call.variant for call in recall_calls) != expected_keys:
        raise ValueError("recall calls must follow the leading candidate allele order")

    lane_ranks: dict[str, int | None] = {}
    for label, ranked in sorted(lanes.items()):
        lane_ranks[label] = next(
            (
                rank
                for rank, candidate in enumerate(ranked.candidates, start=1)
                if _candidate_identity(candidate) == identity
            ),
            None,
        )

    allele_map = {item.variant: item for item in read_audit.alleles}
    allele_reads = tuple(allele_map[key] for key in expected_keys if key in allele_map)
    key_set = frozenset(expected_keys)
    pair_read = next(
        (pair for pair in read_audit.pairs if frozenset((pair.left, pair.right)) == key_set),
        None,
    )
    clinvar_matches = tuple(match for match in clinvar.matches if match.variant in key_set)
    return LeadingCandidateAdjudication(
        leading_candidate=leading,
        lane_ranks=lane_ranks,
        recall_calls=recall_calls,
        recall_status=_recall_status(recall_calls),
        recall_pair_phase=_recall_pair_phase(recall_calls),
        allele_read_evidence=allele_reads,
        pair_read_evidence=pair_read,
        clinvar_release=clinvar.release,
        clinvar_matches=clinvar_matches,
    )


def _escape(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _mane(candidate: CandidateHypothesis, key: VariantKey) -> tuple[str, str, str, str]:
    variant = next(item for item in candidate.variants if item.key == key)
    transcript = next(
        (item for item in variant.annotation.transcripts if item.mane_select),
        next(iter(variant.annotation.transcripts), None),
    )
    if transcript is None:
        return "unavailable", "unavailable", "unavailable", "unavailable"
    return (
        f"{transcript.mane_select or 'unavailable'} / {transcript.transcript or 'unavailable'}",
        transcript.exon or transcript.intron or "unavailable",
        transcript.hgvsc or "unavailable",
        transcript.hgvsp or "unavailable",
    )


def _genotype_summary(call: VcfCallObservation) -> str:
    if not call.present:
        if call.sample is not None:
            return (
                f"exact record present, but genotype {call.gt or 'n/a'} does not contain "
                "the candidate ALT"
            )
        return "not recalled exactly"
    ab = f"{call.allele_balance:.3f}" if call.allele_balance is not None else "n/a"
    return (
        f"{call.gt or 'n/a'}; DP={call.dp}; AD={call.ad_ref},{call.ad_alt}; "
        f"AB={ab}; GQ={call.gq}; QUAL={call.qual}; FILTER={call.filter}"
    )


def _read_summary(read: AlleleReadEvidence | None) -> str:
    if read is None:
        return "not audited"
    ab = f"{read.allele_balance:.3f}" if read.allele_balance is not None else "n/a"
    return (
        f"ref={read.ref_reads}; alt={read.alt_reads}; other={read.other_reads}; "
        f"AB={ab}; alt F/R={read.alt_forward_reads}/{read.alt_reverse_reads}; "
        f"proper-pair={read.alt_proper_pair_reads}; soft-clipped={read.alt_soft_clipped_reads}; "
        f"mean-MQ={read.mean_alt_mapping_quality}; mean-BQ={read.mean_alt_base_quality}; "
        f"mean-end-distance={read.mean_alt_distance_from_read_end}"
    )


def _clinvar_summary(matches: tuple[ClinVarMatch, ...]) -> str:
    if not matches:
        return "no exact match in queried release"
    return "; ".join(
        f"{match.clinical_significance} ({match.review_status}; "
        f"submitters={match.number_submitters if match.number_submitters is not None else 'n/a'}; "
        f"last={match.last_evaluated or 'n/a'}; VariationID={match.variation_id}; "
        f"RCV={','.join(match.rcv_accessions) or 'n/a'}; "
        f"conditions={','.join(match.phenotypes) or 'n/a'})"
        for match in matches
    )


def render_leading_adjudication(adjudication: LeadingCandidateAdjudication) -> str:
    candidate = adjudication.leading_candidate
    read_by_key = {item.variant: item for item in adjudication.allele_read_evidence}
    call_by_key = {item.variant: item for item in adjudication.recall_calls}
    clinvar_by_key = {
        key: tuple(match for match in adjudication.clinvar_matches if match.variant == key)
        for key in (variant.key for variant in candidate.variants)
    }
    alleles = " + ".join(f"`{variant.key.label}`" for variant in candidate.variants)
    lines = [
        "# Private leading-candidate adjudication",
        "",
        "> Research prioritization only. This report neither diagnoses the participant "
        "nor creates a submission.",
        "",
        "## Decision snapshot",
        "",
        f"- Gene / inheritance: `{candidate.gene}` / `{candidate.inheritance}`",
        f"- Exact hypothesis: {alleles}",
        f"- Focused score / ordinal EPCR: `{candidate.score:.6f}` / `{candidate.epcr:.9f}`",
        f"- Exact recall status: `{adjudication.recall_status}`",
        f"- WhatsHap/recall phase: `{adjudication.recall_pair_phase}`",
        "- Direct fragment phase: `"
        + (
            adjudication.pair_read_evidence.interpretation
            if adjudication.pair_read_evidence
            else "not_audited"
        )
        + "`",
        "- Classification policy: no automatic ACMG/AMP classification; evidence "
        "categories remain separate for manual review",
        "",
        "## Cross-lane rank",
        "",
        "| Lane | Exact-pair rank |",
        "|---|---:|",
    ]
    for lane, rank in adjudication.lane_ranks.items():
        lines.append(f"| {_escape(lane)} | {f'#{rank}' if rank is not None else 'not ranked'} |")

    lines.extend(
        [
            "",
            "## Exact variant evidence matrix",
            "",
            "| Variant | MANE RefSeq / Ensembl | Exon/intron | HGVSc | HGVSp | Consequence | "
            "Supplied call | Recall call | Direct BAM | Population AF | Exact ClinVar | "
            "Predictors |",
            "|---|---|---|---|---|---|---|---|---|---:|---|---|",
        ]
    )
    for variant in candidate.variants:
        key = variant.key
        transcript, exon_or_intron, hgvsc, hgvsp = _mane(candidate, key)
        supplied = variant.genotype
        supplied_ab = (
            f"{supplied.allele_balance:.3f}" if supplied.allele_balance is not None else "n/a"
        )
        supplied_summary = (
            f"{supplied.gt}; DP={supplied.dp}; AD={supplied.ad_ref},{supplied.ad_alt}; "
            f"AB={supplied_ab}; GQ={supplied.gq}; QUAL={supplied.qual}; "
            f"FILTER={supplied.filter}"
        )
        annotation = variant.annotation
        population = (
            f"{annotation.max_population_af:.6g}"
            if annotation.max_population_af is not None
            else "unavailable"
        )
        am = annotation.alphamissense if annotation.alphamissense is not None else "n/a"
        sift = (
            annotation.sift_deleteriousness
            if annotation.sift_deleteriousness is not None
            else "n/a"
        )
        polyphen = (
            annotation.polyphen_damagingness
            if annotation.polyphen_damagingness is not None
            else "n/a"
        )
        predictors = f"AM={am}; SIFT-del={sift}; PolyPhen-dmg={polyphen}"
        values = (
            f"`{key.label}`",
            f"`{transcript}`",
            f"`{exon_or_intron}`",
            f"`{hgvsc}`",
            f"`{hgvsp}`",
            annotation.consequence,
            supplied_summary,
            _genotype_summary(call_by_key[key]),
            _read_summary(read_by_key.get(key)),
            population,
            _clinvar_summary(clinvar_by_key[key]),
            predictors,
        )
        lines.append("| " + " | ".join(_escape(value) for value in values) + " |")

    phenotype_scores = {variant.annotation.phenotype_gene_score for variant in candidate.variants}
    mechanism_scores = {
        variant.annotation.disease_mechanism_match for variant in candidate.variants
    }
    phenotype_score_text = ", ".join(str(value) for value in sorted(phenotype_scores))
    mechanism_score_text = ", ".join(str(value) for value in sorted(mechanism_scores))
    lines.extend(
        [
            "",
            "## Pair-level evidence",
            "",
            f"- Focused-lane phase: `{candidate.phase_status.value}`",
            f"- Phenotype-gene score(s): `{phenotype_score_text}`",
            f"- Disease-mechanism score(s): `{mechanism_score_text}`",
            f"- ClinVar exact-match release: `{adjudication.clinvar_release}`",
        ]
    )
    if adjudication.pair_read_evidence is not None:
        pair = adjudication.pair_read_evidence
        lines.append(
            "- Shared-fragment counts: "
            f"alt/alt={pair.alt_alt_fragments}, alt/ref={pair.alt_ref_fragments}, "
            f"ref/alt={pair.ref_alt_fragments}, ref/ref={pair.ref_ref_fragments}, "
            f"other={pair.other_fragments}"
        )

    lines.extend(
        [
            "",
            "## Manual interpretation worksheet",
            "",
            "| Evidence category | What this report establishes | What remains unproven |",
            "|---|---|---|",
            "| Technical allele support | Exact supplied, recall, and BAM observations are "
            "shown independently. | Caller concordance uses the same reads and reference; "
            "it is not biological replication. |",
            "| Population evidence | Exact available maximum AF values are shown; missing "
            "remains unavailable. | Missing frequency is not proof that an allele is absent "
            "from populations. |",
            "| Clinical assertions | Only exact GRCh38 ClinVar allele matches are shown with "
            "review status. | Assertions for another nucleotide or merely the same protein "
            "change are not transferred. |",
            "| Computational evidence | AlphaMissense, SIFT, and PolyPhen values are exposed "
            "separately. | Correlated predictors do not prove pathogenicity and must not be "
            "counted as independent votes. |",
            "| Loss-of-function evidence | Transcript consequences and MANE notation are "
            "explicit. | Formal PVS1 strength requires manual transcript, NMD, mechanism, "
            "and rescue review. |",
            "| Phase / inheritance | Supplied, WhatsHap, and spanning-fragment phase evidence "
            "are separated. | Without a shared informative block or family segregation, "
            "trans configuration remains unresolved. |",
            "| Phenotype / mechanism | Versioned phenotype-gene and mechanism scores are "
            "shown. | Phenotype fit is supportive, not variant-specific functional evidence. |",
            "",
            "## Recorded cautions",
            "",
        ]
    )
    lines.extend(f"- {caution}" for caution in candidate.cautions)
    if not candidate.cautions:
        lines.append("- None recorded by the ranking policy.")
    return "\n".join(lines).rstrip() + "\n"


def write_leading_adjudication(
    adjudication: LeadingCandidateAdjudication, *, json_output: Path, review_output: Path
) -> None:
    json_output.parent.mkdir(parents=True, exist_ok=True)
    review_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(adjudication.model_dump_json(indent=2) + "\n", encoding="utf-8")
    review_output.write_text(render_leading_adjudication(adjudication), encoding="utf-8")
