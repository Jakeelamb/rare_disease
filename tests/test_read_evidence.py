from pathlib import Path

import pysam

from mva_hackathon.models import GenotypeEvidence, VariantAnnotation, VariantEvidence, VariantKey
from mva_hackathon.policy import load_policy
from mva_hackathon.ranking import rank_case
from mva_hackathon.read_evidence import (
    PairReadEvidence,
    inspect_ranked_bam,
    render_candidate_read_audit,
)


def _ranked_pair():
    variants = [
        VariantEvidence(
            key=VariantKey(chrom="1", pos=position, ref=ref, alt=alt),
            genotype=GenotypeEvidence(gt="0/1", dp=40, ad_ref=20, ad_alt=20, gq=99),
            annotation=VariantAnnotation(
                gene="SYNTHETIC",
                consequence="frameshift_variant",
                max_population_af=0.00001,
            ),
        )
        for position, ref, alt in ((100, "A", "G"), (110, "C", "T"))
    ]
    return rank_case(variants, load_policy(Path("config/scoring_policy.yaml")))


def _write_bam(path: Path) -> None:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": 1000}]}
    patterns = (
        ("alt_ref_1", "G", "C"),
        ("alt_ref_2", "G", "C"),
        ("ref_alt_1", "A", "T"),
        ("ref_alt_2", "A", "T"),
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, first, second in patterns:
            sequence = list("A" * 100)
            sequence[49] = first
            sequence[59] = second
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "".join(sequence)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 50
            read.mapping_quality = 60
            read.cigar = ((0, 100),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)
    pysam.index(str(path))


def test_candidate_read_audit_counts_alleles_and_read_backed_trans(tmp_path: Path) -> None:
    bam = tmp_path / "synthetic.bam"
    _write_bam(bam)

    audit = inspect_ranked_bam(_ranked_pair(), bam)

    assert len(audit.alleles) == 2
    assert [(item.ref_reads, item.alt_reads) for item in audit.alleles] == [(2, 2), (2, 2)]
    assert all(item.allele_balance == 0.5 for item in audit.alleles)
    assert all(item.alt_proper_pair_reads == 0 for item in audit.alleles)
    assert all(item.alt_soft_clipped_reads == 0 for item in audit.alleles)
    assert [item.mean_alt_distance_from_read_end for item in audit.alleles] == [49.0, 40.0]
    assert len(audit.pairs) == 1
    assert audit.pairs[0].rank == 1
    assert audit.pairs[0].gene == "SYNTHETIC"
    assert audit.pairs[0].shared_fragments == 4
    assert audit.pairs[0].min_phase_support_fragments == 2
    assert audit.pairs[0].alt_ref_fragments == 2
    assert audit.pairs[0].ref_alt_fragments == 2
    assert audit.pairs[0].interpretation == "read_backed_trans"
    assert "read_backed_trans" in render_candidate_read_audit(audit)
    assert "| 1 | SYNTHETIC |" in render_candidate_read_audit(audit)
    assert "Mean end distance" in render_candidate_read_audit(audit)


def test_ref_only_spanning_fragment_is_not_called_conflicting_phase() -> None:
    left = VariantKey(chrom="1", pos=100, ref="A", alt="G")
    right = VariantKey(chrom="1", pos=110, ref="C", alt="T")

    evidence = PairReadEvidence(
        left=left,
        right=right,
        shared_fragments=1,
        alt_alt_fragments=0,
        alt_ref_fragments=0,
        ref_alt_fragments=0,
        ref_ref_fragments=1,
        other_fragments=0,
    )

    assert evidence.interpretation == "no_informative_spanning_fragments"


def test_pair_phase_requires_exact_alt_support_for_both_alleles() -> None:
    left = VariantKey(chrom="1", pos=100, ref="A", alt="AG")
    right = VariantKey(chrom="1", pos=110, ref="C", alt="T")

    evidence = PairReadEvidence(
        left=left,
        right=right,
        left_alt_observed=False,
        right_alt_observed=True,
        shared_fragments=4,
        alt_alt_fragments=0,
        alt_ref_fragments=0,
        ref_alt_fragments=4,
        ref_ref_fragments=0,
        other_fragments=0,
    )

    assert evidence.interpretation == "insufficient_exact_alt_support"


def test_single_concordant_fragment_is_support_not_a_phase_call() -> None:
    left = VariantKey(chrom="1", pos=100, ref="A", alt="G")
    right = VariantKey(chrom="1", pos=110, ref="C", alt="T")

    evidence = PairReadEvidence(
        left=left,
        right=right,
        shared_fragments=1,
        alt_alt_fragments=1,
        alt_ref_fragments=0,
        ref_alt_fragments=0,
        ref_ref_fragments=0,
        other_fragments=0,
    )

    assert evidence.interpretation == "insufficient_phase_support"
