from mva_hackathon.models import VariantKey
from mva_hackathon.scorer import ScoredRow, score_rows


def key(chrom: str, pos: int) -> VariantKey:
    return VariantKey(chrom=chrom, pos=pos, ref="A", alt="G")


def test_exact_pair_at_rank_one_gets_full_score() -> None:
    truth = frozenset((key("1", 100), key("1", 200)))
    result = score_rows([ScoredRow(truth, 0.9, 99)], truth)

    assert result.full_match_rank == 1
    assert result.rank_points == 100
    assert result.f_max == 1.0


def test_split_true_alleles_get_partial_rank_but_full_fmax() -> None:
    left, right = key("1", 100), key("1", 200)
    truth = frozenset((left, right))
    rows = [
        ScoredRow(frozenset((left,)), 0.9, 1),
        ScoredRow(frozenset((right,)), 0.8, 2),
    ]

    result = score_rows(rows, truth)

    assert result.full_match_rank is None
    assert result.partial_match_rank == 1
    assert result.rank_points == 50
    assert result.f_max == 1.0
    assert result.f_max_threshold == 0.8


def test_false_variant_above_pair_reduces_rank_and_fmax() -> None:
    truth = frozenset((key("1", 100), key("1", 200)))
    rows = [
        ScoredRow(frozenset((key("2", 300),)), 0.9, 1),
        ScoredRow(truth, 0.8, 2),
    ]

    result = score_rows(rows, truth)

    assert result.full_match_rank == 2
    assert result.rank_points == 50
    assert result.f_max == 0.8
