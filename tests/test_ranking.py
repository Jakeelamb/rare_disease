from pathlib import Path

import pytest

from mva_hackathon.models import (
    GenotypeEvidence,
    PhaseStatus,
    VariantAnnotation,
    VariantEvidence,
    VariantKey,
)
from mva_hackathon.policy import load_policy
from mva_hackathon.ranking import rank_case

POLICY = load_policy(Path("config/scoring_policy.yaml"))


def variant(
    pos: int,
    *,
    gene: str = "GENE1",
    gt: str = "0/1",
    phased_gt: str | None = None,
    phase_set: str | None = None,
    consequence: str = "frameshift_variant",
    af: float | None = 0.00001,
    phenotype_score: float = 0.9,
    cadd_phred: float | None = 30,
    dp: int | None = 40,
    gq: float | None = 99,
    clinvar_significance: str | None = None,
) -> VariantEvidence:
    return VariantEvidence(
        key=VariantKey(chrom="15", pos=pos, ref="A", alt="T"),
        genotype=GenotypeEvidence(
            gt=gt,
            dp=dp,
            ad_ref=20,
            ad_alt=20,
            gq=gq,
            phase_set=phase_set,
            phased_gt=phased_gt,
            phase_method="read_backed" if phased_gt else None,
        ),
        annotation=VariantAnnotation(
            gene=gene,
            consequence=consequence,
            max_population_af=af,
            phenotype_gene_score=phenotype_score,
            disease_mechanism_match=0.9,
            cadd_phred=cadd_phred,
            clinvar_significance=clinvar_significance,
        ),
    )


def test_builds_exact_same_gene_compound_heterozygous_pair() -> None:
    result = rank_case([variant(100), variant(200)], POLICY)

    assert len(result.candidates) == 1
    top = result.candidates[0]
    assert top.inheritance == "compound_heterozygous"
    assert {item.key.pos for item in top.variants} == {100, 200}
    assert top.phase_status == PhaseStatus.UNRESOLVED
    assert top.cautions == ("trans phase is unresolved",)
    assert [item.term for item in top.contributions][-2:] == ["pair", "phase"]


def test_read_backed_cis_pair_is_excluded() -> None:
    left = variant(100, phase_set="block1", phased_gt="0|1")
    right = variant(200, phase_set="block1", phased_gt="0|1")

    result = rank_case([left, right], POLICY)

    assert result.candidates == ()


def test_read_backed_trans_pair_gets_phase_support() -> None:
    left = variant(100, phase_set="block1", phased_gt="0|1")
    right = variant(200, phase_set="block1", phased_gt="1|0")

    top = rank_case([left, right], POLICY).candidates[0]

    assert top.phase_status == PhaseStatus.READ_BACKED_TRANS
    assert next(item for item in top.contributions if item.term == "phase").raw_value == 0.9


def test_phased_genotypes_without_evidence_method_remain_unresolved() -> None:
    left = variant(100, phase_set="block1", phased_gt="0|1").model_copy(
        update={
            "genotype": GenotypeEvidence(
                gt="0/1",
                dp=40,
                ad_ref=20,
                ad_alt=20,
                gq=99,
                phase_set="block1",
                phased_gt="0|1",
            )
        }
    )
    right = variant(200, phase_set="block1", phased_gt="1|0").model_copy(
        update={
            "genotype": GenotypeEvidence(
                gt="0/1",
                dp=40,
                ad_ref=20,
                ad_alt=20,
                gq=99,
                phase_set="block1",
                phased_gt="1|0",
            )
        }
    )

    assert rank_case([left, right], POLICY).candidates[0].phase_status == PhaseStatus.UNRESOLVED


def test_low_quality_or_common_alleles_are_excluded() -> None:
    common = variant(100, af=0.2)
    low_balance = variant(200).model_copy(
        update={"genotype": GenotypeEvidence(gt="0/1", dp=40, ad_ref=39, ad_alt=1, gq=99)}
    )

    result = rank_case([common, low_balance], POLICY)

    assert result.candidates == ()
    assert result.excluded_variant_count == 2


def test_homozygous_alternative_is_preserved() -> None:
    result = rank_case([variant(100, gt="1/1")], POLICY)

    assert result.candidates[0].inheritance == "recessive_homozygous"


def test_epcr_is_strictly_descending_and_bounded() -> None:
    values = [
        variant(100),
        variant(200),
        variant(300, consequence="missense_variant"),
    ]

    probabilities = [item.epcr for item in rank_case(values, POLICY).candidates]

    assert probabilities == sorted(probabilities, reverse=True)
    assert len(probabilities) == len(set(probabilities))
    assert all(0 < probability <= 1 for probability in probabilities)


def test_missing_evidence_never_receives_positive_credit() -> None:
    missing = variant(100, af=None, cadd_phred=None, dp=None, gq=None)
    complete = variant(200)

    top = rank_case([missing, complete], POLICY).candidates[0]
    contributions = {item.term: item for item in top.contributions}

    assert contributions["rarity"].raw_value == 0.0
    assert contributions["pathogenicity"].raw_value == 0.0
    assert "missing evidence" in " ".join(top.cautions)


@pytest.mark.parametrize(
    ("label", "expected"),
    (("Likely pathogenic", 0.85), ("Likely benign", -0.70)),
)
def test_clinvar_labels_match_exact_policy_categories(label: str, expected: float) -> None:
    top = rank_case([variant(100, gt="1/1", clinvar_significance=label)], POLICY).candidates[0]
    clinical = next(item for item in top.contributions if item.term == "clinical")

    assert clinical.raw_value == expected


def test_loss_of_function_alleles_do_not_request_missense_predictors() -> None:
    left = variant(100, consequence="stop_gained", cadd_phred=None)
    right = variant(200, consequence="frameshift_variant", cadd_phred=None)

    top = rank_case([left, right], POLICY).candidates[0]

    assert "pathogenicity predictors" not in " ".join(top.cautions)


def test_weak_passenger_allele_cannot_form_a_pair() -> None:
    strong = variant(100)
    passenger = variant(
        200,
        consequence="synonymous_variant",
        cadd_phred=None,
        phenotype_score=0.9,
    )

    result = rank_case([strong, passenger], POLICY)

    assert result.candidates == ()
    assert result.low_support_variant_count == 1


def test_pair_terms_use_weakest_allele_not_mean() -> None:
    strong = variant(100, consequence="frameshift_variant")
    weaker = variant(200, consequence="missense_variant")

    top = rank_case([strong, weaker], POLICY).candidates[0]
    consequence = next(item for item in top.contributions if item.term == "consequence")

    assert consequence.raw_value == POLICY.consequence_scores["missense_variant"]
    assert "weakest allele" in consequence.rationale


def test_variant_limit_fails_loudly_instead_of_silent_truncation() -> None:
    thresholds = POLICY.thresholds.model_copy(update={"max_variants_per_gene": 2})
    bounded_policy = POLICY.model_copy(update={"thresholds": thresholds})

    with pytest.raises(ValueError, match=r"eligible variants.*configured safety limit"):
        rank_case([variant(100), variant(200), variant(300)], bounded_policy)


def test_variant_key_normalizes_only_representation_not_alleles() -> None:
    key = VariantKey(chrom="chr15", pos=123, ref="a", alt="t")

    assert key.chrom == "15"
    assert key.challenge_chrom == "chr15"
    assert key.label == "chr15:123:A>T"


def test_policy_rejects_missing_required_weight(tmp_path: Path) -> None:
    broken = tmp_path / "broken.yaml"
    text = Path("config/scoring_policy.yaml").read_text().replace("  phase: 10.0\n", "")
    broken.write_text(text)

    with pytest.raises(ValueError, match="missing weights"):
        load_policy(broken)
