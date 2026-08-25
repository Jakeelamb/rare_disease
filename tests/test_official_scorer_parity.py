import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

from mva_hackathon.models import VariantKey
from mva_hackathon.scorer import ScoredRow, score_rows


def _load_official() -> ModuleType:
    path = Path("challenge_space/evaluation.py")
    if not path.exists():
        pytest.skip("pinned official Space checkout is not present")
    spec = importlib.util.spec_from_file_location("official_evaluation", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_local_scorer_matches_pinned_official_evaluator() -> None:
    official = _load_official()
    left = ("chr1", 100, "A", "G")
    right = ("chr1", 200, "C", "T")
    distractor = ("chr2", 300, "G", "A")
    official_rows = [
        official.SubmissionRow(frozenset((distractor,)), 0.95, 1),
        official.SubmissionRow(frozenset((left, right)), 0.80, 2),
    ]
    expected = official.score_proband("PROBAND01", official_rows, frozenset((left, right)))

    local_rows = [
        ScoredRow(
            frozenset((VariantKey(chrom="2", pos=300, ref="G", alt="A"),)),
            0.95,
            1,
        ),
        ScoredRow(
            frozenset(
                (
                    VariantKey(chrom="1", pos=100, ref="A", alt="G"),
                    VariantKey(chrom="1", pos=200, ref="C", alt="T"),
                )
            ),
            0.80,
            2,
        ),
    ]
    observed = score_rows(local_rows, frozenset(local_rows[1].variants))

    assert observed.full_match_rank == expected.full_match_rank
    assert observed.partial_match_rank == expected.partial_match_rank
    assert observed.rank_points == expected.rank_points
    assert observed.f_max == expected.f_max
    assert observed.f_max_threshold == expected.f_max_threshold
    assert observed.n_predictions_at_f_max == expected.n_predictions_at_f_max
