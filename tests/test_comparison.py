from mva_hackathon.comparison import compare_ranked_cases
from mva_hackathon.models import (
    CandidateHypothesis,
    GenotypeEvidence,
    RankedCase,
    VariantAnnotation,
    VariantEvidence,
    VariantKey,
)


def candidate(gene: str, positions: tuple[int, int], score: float) -> CandidateHypothesis:
    variants = tuple(
        VariantEvidence(
            key=VariantKey(chrom="1", pos=position, ref="A", alt="G"),
            genotype=GenotypeEvidence(gt="0/1", dp=30, ad_ref=15, ad_alt=15, gq=99),
            annotation=VariantAnnotation(gene=gene, consequence="frameshift_variant"),
        )
        for position in positions
    )
    return CandidateHypothesis(
        variants=variants,
        gene=gene,
        inheritance="compound_heterozygous",
        score=score,
        epcr=0.5,
        contributions=(),
    )


def ranked(*candidates: CandidateHypothesis) -> RankedCase:
    return RankedCase(
        proband_id="PROBAND01",
        policy_version="test",
        candidates=candidates,
    )


def test_comparison_uses_union_with_lane_provenance_not_majority_vote() -> None:
    shared = candidate("SHARED", (100, 200), 10)
    rescue = candidate("RESCUE", (300, 400), 9)

    comparison = compare_ranked_cases(
        {"supplied": ranked(shared), "deepvariant": ranked(shared, rescue)}
    )

    assert comparison.hypothesis_count == 2
    rescue_row = next(item for item in comparison.hypotheses if item.gene == "RESCUE")
    assert set(rescue_row.lanes) == {"deepvariant"}
    assert rescue_row.review_required
    assert comparison.consensus_count == 1
