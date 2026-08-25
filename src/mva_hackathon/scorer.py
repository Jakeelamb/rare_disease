"""Revision-pinned local mirror of the public Track 1 scoring semantics."""

from __future__ import annotations

from dataclasses import dataclass

from .models import VariantKey


@dataclass(frozen=True)
class ScoredRow:
    variants: frozenset[VariantKey]
    epcr: float
    rank: int


@dataclass(frozen=True)
class ScoreResult:
    full_match_rank: int | None
    partial_match_rank: int | None
    rank_points: float
    f_max: float
    f_max_threshold: float | None
    n_predictions_at_f_max: int


def _rank_points(rank: int) -> int:
    for ceiling, points in ((1, 100), (3, 50), (5, 25), (10, 10)):
        if rank <= ceiling:
            return points
    return 0


def score_rows(rows: list[ScoredRow], true_variants: frozenset[VariantKey]) -> ScoreResult:
    """Score synthetic/review keys exactly as the pinned public evaluator does."""
    ordered = sorted(enumerate(rows), key=lambda item: (-item[1].epcr, item[0]))
    ranked = [
        ScoredRow(variants=row.variants, epcr=row.epcr, rank=index + 1)
        for index, (_, row) in enumerate(ordered)
    ]

    full_match_rank = next((row.rank for row in ranked if row.variants == true_variants), None)
    partial_match_rank = None
    if len(true_variants) == 2 and full_match_rank is None:
        partial_match_rank = next(
            (row.rank for row in ranked if row.variants & true_variants), None
        )

    if full_match_rank is not None:
        rank_points = float(_rank_points(full_match_rank))
    elif partial_match_rank is not None:
        rank_points = 0.5 * _rank_points(partial_match_rank)
    else:
        rank_points = 0.0

    best_f = 0.0
    best_threshold = None
    best_rows = 0
    for threshold in sorted({row.epcr for row in ranked}, reverse=True):
        accepted = [row for row in ranked if row.epcr >= threshold]
        predicted = set().union(*(row.variants for row in accepted))
        true_positive = len(predicted & true_variants)
        false_positive = len(predicted - true_variants)
        false_negative = len(true_variants - predicted)
        precision = (
            true_positive / (true_positive + false_positive)
            if true_positive + false_positive
            else 0.0
        )
        recall = (
            true_positive / (true_positive + false_negative)
            if true_positive + false_negative
            else 0.0
        )
        f_score = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        if f_score > best_f:
            best_f = f_score
            best_threshold = threshold
            best_rows = len(accepted)

    return ScoreResult(
        full_match_rank=full_match_rank,
        partial_match_rank=partial_match_rank,
        rank_points=rank_points,
        f_max=best_f,
        f_max_threshold=best_threshold,
        n_predictions_at_f_max=best_rows,
    )
