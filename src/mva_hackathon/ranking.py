"""Inheritance-aware, fully decomposed candidate ranking."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from itertools import combinations

from .models import (
    CandidateHypothesis,
    PhaseStatus,
    RankedCase,
    ScoreContribution,
    VariantEvidence,
)
from .policy import ScoringPolicy


def _clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(value, high))


def _contribution(
    term: str, raw: float, policy: ScoringPolicy, rationale: str
) -> ScoreContribution:
    weight = policy.weights[term]
    return ScoreContribution(
        term=term,
        raw_value=round(raw, 6),
        weight=weight,
        points=round(raw * weight, 6),
        rationale=rationale,
    )


def _quality(variant: VariantEvidence, policy: ScoringPolicy) -> ScoreContribution:
    gt = variant.genotype
    parts = [1.0 if gt.filter in {"PASS", "."} else 0.0]
    parts.append(0.0 if gt.dp is None else _clamp(gt.dp / policy.thresholds.min_depth))
    parts.append(0.0 if gt.gq is None else _clamp(gt.gq / policy.thresholds.min_gq))
    balance = gt.allele_balance
    if gt.is_heterozygous:
        if balance is None:
            parts.append(0.0)
        else:
            center_distance = abs(balance - 0.5) / 0.5
            parts.append(_clamp(1.0 - center_distance))
    else:
        parts.append(1.0 if gt.is_homozygous_alt else 0.25)
    raw = sum(parts) / len(parts)
    return _contribution(
        "quality",
        raw,
        policy,
        f"FILTER={gt.filter}; DP={gt.dp}; GQ={gt.gq}; AB={balance}",
    )


def _rarity(variant: VariantEvidence, policy: ScoringPolicy) -> ScoreContribution:
    af = variant.annotation.max_population_af
    if af is None:
        raw, rationale = 0.0, "population frequency unavailable; no positive credit"
    elif af <= 1e-5:
        raw, rationale = 1.0, f"max population AF {af:g} <= 1e-5"
    elif af <= 1e-4:
        raw, rationale = 0.9, f"max population AF {af:g} <= 1e-4"
    elif af <= 1e-3:
        raw, rationale = 0.7, f"max population AF {af:g} <= 1e-3"
    elif af <= policy.thresholds.max_recessive_af:
        raw, rationale = 0.2, f"max population AF {af:g} within recessive ceiling"
    else:
        raw, rationale = 0.0, f"max population AF {af:g} exceeds recessive ceiling"
    return _contribution("rarity", raw, policy, rationale)


def _consequence(variant: VariantEvidence, policy: ScoringPolicy) -> ScoreContribution:
    consequences = variant.annotation.consequence.lower().split("&")
    scores = [policy.consequence_scores.get(item, 0.0) for item in consequences]
    raw = max(scores, default=0.0)
    return _contribution(
        "consequence",
        raw,
        policy,
        f"most severe annotated consequence: {variant.annotation.consequence}",
    )


def _clinical(variant: VariantEvidence, policy: ScoringPolicy) -> ScoreContribution:
    significance = (variant.annotation.clinvar_significance or "unknown").lower()
    normalized_significance = re.sub(r"[^a-z0-9]+", "_", significance).strip("_")
    matched = next(
        (
            score
            for label, score in sorted(
                policy.clinvar_scores.items(), key=lambda item: len(item[0]), reverse=True
            )
            if label in normalized_significance
        ),
        policy.clinvar_scores.get("unknown", 0.0),
    )
    return _contribution(
        "clinical",
        matched,
        policy,
        f"ClinVar={variant.annotation.clinvar_significance or 'unavailable'}; "
        f"review={variant.annotation.clinvar_review_status or 'unavailable'}",
    )


def _pathogenicity(variant: VariantEvidence, policy: ScoringPolicy) -> ScoreContribution:
    ann = variant.annotation
    values: list[tuple[str, float]] = []
    if ann.cadd_phred is not None:
        values.append(("CADD", _clamp((ann.cadd_phred - 10.0) / 20.0)))
    if ann.revel is not None:
        values.append(("REVEL", ann.revel))
    if ann.spliceai is not None:
        values.append(("SpliceAI", ann.spliceai))
    if ann.alphamissense is not None:
        values.append(("AlphaMissense", ann.alphamissense))
    native_missense = [
        value
        for value in (ann.sift_deleteriousness, ann.polyphen_damagingness)
        if value is not None
    ]
    if native_missense:
        values.append(("VEP missense consensus", sum(native_missense) / len(native_missense)))
    if not values:
        return _contribution(
            "pathogenicity",
            0.0,
            policy,
            "computational predictions unavailable; no positive credit",
        )
    raw = sum(value for _, value in values) / len(values)
    detail = ", ".join(f"{name}={value:.3f}" for name, value in values)
    return _contribution("pathogenicity", raw, policy, detail)


def _variant_contributions(
    variant: VariantEvidence, policy: ScoringPolicy
) -> tuple[ScoreContribution, ...]:
    ann = variant.annotation
    return (
        _quality(variant, policy),
        _rarity(variant, policy),
        _consequence(variant, policy),
        _clinical(variant, policy),
        _pathogenicity(variant, policy),
        _contribution(
            "phenotype_gene",
            ann.phenotype_gene_score,
            policy,
            f"versioned phenotype-gene score for {ann.gene}",
        ),
        _contribution(
            "mechanism",
            ann.disease_mechanism_match,
            policy,
            f"curated disease-mechanism match for {ann.gene}",
        ),
    )


def _infer_phase(left: VariantEvidence, right: VariantEvidence) -> PhaseStatus:
    lgt, rgt = left.genotype, right.genotype
    if not lgt.phase_method or lgt.phase_method != rgt.phase_method:
        return PhaseStatus.UNRESOLVED
    if not lgt.phase_set or lgt.phase_set != rgt.phase_set:
        return PhaseStatus.UNRESOLVED
    if not lgt.phased_gt or not rgt.phased_gt or "|" not in lgt.phased_gt + rgt.phased_gt:
        return PhaseStatus.UNRESOLVED

    def alt_haplotype(phased_gt: str) -> int | None:
        alleles = phased_gt.split("|")
        return alleles.index("1") if len(alleles) == 2 and alleles.count("1") == 1 else None

    left_hap = alt_haplotype(lgt.phased_gt)
    right_hap = alt_haplotype(rgt.phased_gt)
    if left_hap is None or right_hap is None:
        return PhaseStatus.UNRESOLVED
    if left_hap == right_hap:
        return PhaseStatus.INCOMPATIBLE_CIS
    return {
        "family": PhaseStatus.FAMILY_BACKED_TRANS,
        "read_backed": PhaseStatus.READ_BACKED_TRANS,
        "statistical": PhaseStatus.STATISTICAL_TRANS,
    }[lgt.phase_method]


def _eligible(variant: VariantEvidence, policy: ScoringPolicy) -> bool:
    gt = variant.genotype
    af = variant.annotation.max_population_af
    balance = gt.allele_balance
    balance_ok = (
        not gt.is_heterozygous
        or balance is None
        or policy.thresholds.min_het_allele_balance
        <= balance
        <= policy.thresholds.max_het_allele_balance
    )
    return (
        gt.filter in {"PASS", "."}
        and (gt.dp is None or gt.dp >= policy.thresholds.min_depth)
        and (gt.gq is None or gt.gq >= policy.thresholds.min_gq)
        and (af is None or af <= policy.thresholds.max_recessive_af)
        and (gt.is_heterozygous or gt.is_homozygous_alt)
        and balance_ok
    )


def _allele_support(contributions: tuple[ScoreContribution, ...]) -> float:
    values = {item.term: item.raw_value for item in contributions}
    return max(
        values["consequence"],
        max(values["clinical"], 0.0),
        values["pathogenicity"],
    )


def _missing_evidence_caution(variant: VariantEvidence) -> str | None:
    missing: list[str] = []
    genotype = variant.genotype
    annotation = variant.annotation
    if genotype.dp is None:
        missing.append("DP")
    if genotype.gq is None:
        missing.append("GQ")
    if genotype.is_heterozygous and genotype.allele_balance is None:
        missing.append("allele balance")
    if annotation.max_population_af is None:
        missing.append("population AF")
    consequences = set(annotation.consequence.lower().split("&"))
    predictor_exempt = bool(
        consequences
        & {
            "transcript_ablation",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "frameshift_variant",
            "start_lost",
        }
    )
    if not predictor_exempt and all(
        value is None
        for value in (
            annotation.cadd_phred,
            annotation.revel,
            annotation.spliceai,
            annotation.alphamissense,
            annotation.sift_deleteriousness,
            annotation.polyphen_damagingness,
        )
    ):
        missing.append("pathogenicity predictors")
    if not missing:
        return None
    return f"{variant.key.label}: missing evidence: {', '.join(missing)}"


def _candidate_sort_key(candidate: CandidateHypothesis) -> tuple[float, str, tuple[str, ...]]:
    labels = tuple(sorted(variant.key.label for variant in candidate.variants))
    return (-candidate.score, candidate.gene, labels)


def _assign_epcr(
    candidates: list[CandidateHypothesis], temperature: float
) -> list[CandidateHypothesis]:
    if not candidates:
        return []
    peak = max(candidate.score for candidate in candidates)
    exponents = [math.exp((candidate.score - peak) / temperature) for candidate in candidates]
    denominator = sum(exponents)
    probabilities = [value / denominator for value in exponents]

    ranked: list[CandidateHypothesis] = []
    previous = 1.0
    for index, (candidate, probability) in enumerate(zip(candidates, probabilities, strict=True)):
        # The evaluator sorts on EPCR. Keep it strictly descending even when scores tie.
        epcr = min(probability, previous - 1e-9) if index else min(probability, 0.999999)
        epcr = max(epcr, 1e-9)
        previous = epcr
        ranked.append(candidate.model_copy(update={"epcr": round(epcr, 9)}))
    return ranked


def rank_case(
    variants: list[VariantEvidence], policy: ScoringPolicy, proband_id: str = "PROBAND01"
) -> RankedCase:
    """Rank one case through the public analysis interface.

    The implementation builds exact same-gene heterozygous pairs, preserves homozygous
    alternatives, exposes every score term, and never treats missing annotation as positive proof.
    """

    eligible = [variant for variant in variants if _eligible(variant, policy)]
    by_gene: dict[str, list[tuple[VariantEvidence, tuple[ScoreContribution, ...], float]]] = (
        defaultdict(list)
    )
    low_support_variant_count = 0
    for variant in eligible:
        contributions = _variant_contributions(variant, policy)
        if _allele_support(contributions) < policy.thresholds.min_allele_support:
            low_support_variant_count += 1
            continue
        score = sum(item.points for item in contributions)
        by_gene[variant.annotation.gene].append((variant, contributions, score))

    candidates: list[CandidateHypothesis] = []
    for gene, records in by_gene.items():
        records.sort(key=lambda item: (-item[2], item[0].key.label))
        if len(records) > policy.thresholds.max_variants_per_gene:
            raise ValueError(
                f"gene {gene} has {len(records)} eligible variants, exceeding configured "
                f"safety limit {policy.thresholds.max_variants_per_gene}; refine or raise the "
                "visible policy bound"
            )

        for variant, contributions, score in records:
            if variant.genotype.is_homozygous_alt:
                missing_caution = _missing_evidence_caution(variant)
                homozygous_cautions = [
                    "single-proband inference; parental segregation is unavailable"
                ]
                if missing_caution:
                    homozygous_cautions.append(missing_caution)
                candidates.append(
                    CandidateHypothesis(
                        variants=(variant,),
                        gene=gene,
                        inheritance="recessive_homozygous",
                        score=score,
                        contributions=contributions,
                        cautions=tuple(homozygous_cautions),
                    )
                )

        heterozygous = [record for record in records if record[0].genotype.is_heterozygous]
        for left, right in combinations(heterozygous, 2):
            left_variant, left_contributions, _ = left
            right_variant, right_contributions, _ = right
            phase = _infer_phase(left_variant, right_variant)
            if phase == PhaseStatus.INCOMPATIBLE_CIS:
                continue
            averaged: list[ScoreContribution] = []
            for left_item, right_item in zip(left_contributions, right_contributions, strict=True):
                gene_level = left_item.term in {"phenotype_gene", "mechanism"}
                raw = (
                    (left_item.raw_value + right_item.raw_value) / 2.0
                    if gene_level
                    else min(left_item.raw_value, right_item.raw_value)
                )
                rule = "gene-level mean" if gene_level else "weakest allele"
                averaged.append(
                    _contribution(
                        left_item.term,
                        raw,
                        policy,
                        f"{rule}: [{left_item.rationale}] + [{right_item.rationale}]",
                    )
                )
            left_support = _allele_support(left_contributions)
            right_support = _allele_support(right_contributions)
            averaged.extend(
                [
                    _contribution(
                        "pair",
                        min(left_support, right_support),
                        policy,
                        "pair support is limited by the weaker eligible allele",
                    ),
                    _contribution(
                        "phase",
                        policy.phase_scores[phase.value],
                        policy,
                        f"phase status: {phase.value}",
                    ),
                ]
            )
            cautions: list[str] = []
            if phase == PhaseStatus.UNRESOLVED:
                cautions.append("trans phase is unresolved")
            for variant in (left_variant, right_variant):
                missing_caution = _missing_evidence_caution(variant)
                if missing_caution:
                    cautions.append(missing_caution)
            candidates.append(
                CandidateHypothesis(
                    variants=(left_variant, right_variant),
                    gene=gene,
                    inheritance="compound_heterozygous",
                    phase_status=phase,
                    score=sum(item.points for item in averaged),
                    contributions=tuple(averaged),
                    cautions=tuple(cautions),
                )
            )

    candidates.sort(key=_candidate_sort_key)
    ranked = _assign_epcr(candidates, policy.epcr_temperature)
    return RankedCase(
        proband_id=proband_id,
        policy_version=policy.version,
        candidates=tuple(ranked),
        excluded_variant_count=len(variants) - len(eligible) + low_support_variant_count,
        low_support_variant_count=low_support_variant_count,
    )
