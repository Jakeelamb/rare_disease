from pathlib import Path

from mva_hackathon.models import (
    GenotypeEvidence,
    TranscriptConsequence,
    VariantAnnotation,
    VariantEvidence,
    VariantKey,
)
from mva_hackathon.policy import load_policy
from mva_hackathon.ranking import rank_case
from mva_hackathon.report import render_ranked_case


def test_report_exposes_pair_score_terms_and_cautions() -> None:
    variants = [
        VariantEvidence(
            key=VariantKey(chrom="1", pos=position, ref="A", alt="G"),
            genotype=GenotypeEvidence(gt="0/1", dp=30, ad_ref=15, ad_alt=15, gq=99),
            annotation=VariantAnnotation(
                gene="SYNTHETIC",
                consequence="frameshift_variant",
                transcripts=(
                    TranscriptConsequence(
                        transcript="ENST_TEST",
                        consequence="frameshift_variant",
                        mane_select="NM_TEST",
                        hgvsc=f"ENST_TEST:c.{position}del",
                        hgvsp="ENSP_TEST:p.Gly1fs",
                    ),
                ),
                max_population_af=0.00001,
                alphamissense=0.9,
            ),
        )
        for position in (100, 200)
    ]
    ranked = rank_case(variants, load_policy(Path("config/scoring_policy.yaml")))

    report = render_ranked_case(ranked)

    assert "`chr1:100:A>G` + `chr1:200:A>G`" in report
    assert "compound_heterozygous" in report
    assert "trans phase is unresolved" in report
    assert "ENST_TEST:c.100del" in report
    assert "ENSP_TEST:p.Gly1fs" in report
    assert "AlphaMissense" in report
    assert "0.9" in report
    assert "| pair |" in report
