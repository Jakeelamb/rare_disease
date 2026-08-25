"""Truth-set metrics and an advisory comparison for additional algorithms."""

from __future__ import annotations

import json
from pathlib import Path

import yaml
from pydantic import BaseModel, Field


class BenchmarkMetrics(BaseModel):
    callset: str
    truth_set: str
    engine: str
    true_positive: int = Field(ge=0)
    false_positive: int = Field(ge=0)
    false_negative: int = Field(ge=0)
    runtime_seconds: float = Field(gt=0)
    peak_rss_mb: float = Field(gt=0)

    @property
    def precision(self) -> float:
        denominator = self.true_positive + self.false_positive
        return self.true_positive / denominator if denominator else 0.0

    @property
    def recall(self) -> float:
        denominator = self.true_positive + self.false_negative
        return self.true_positive / denominator if denominator else 0.0

    @property
    def f1(self) -> float:
        return (
            2 * self.precision * self.recall / (self.precision + self.recall)
            if self.precision + self.recall
            else 0.0
        )


class ComparisonPolicy(BaseModel):
    min_rescued_true_positives: int = Field(default=1, ge=1)
    max_added_false_positives_per_rescue: float = Field(default=1.0, ge=0)
    max_total_runtime_ratio: float = Field(default=3.0, ge=1)


class ToolComparison(BaseModel):
    baseline: str
    candidate: str
    truth_set: str
    meets_policy: bool
    rescued_true_positives: int
    added_false_positives: int
    runtime_ratio: float
    reasons: tuple[str, ...]


def compare_tool_metrics(
    *,
    baseline: BenchmarkMetrics,
    candidate: BenchmarkMetrics,
    union: BenchmarkMetrics,
    policy: ComparisonPolicy,
) -> ToolComparison:
    reasons: list[str] = []
    truth_sets = {baseline.truth_set, candidate.truth_set, union.truth_set}
    engines = {baseline.engine, candidate.engine, union.engine}
    if len(truth_sets) != 1:
        reasons.append("benchmark truth sets differ")
    if len(engines) != 1:
        reasons.append("benchmark engines differ")

    rescued = max(0, union.true_positive - baseline.true_positive)
    added_false_positives = max(0, union.false_positive - baseline.false_positive)
    runtime_ratio = union.runtime_seconds / baseline.runtime_seconds
    if rescued < policy.min_rescued_true_positives:
        reasons.append(
            f"no incremental true-positive rescue meeting minimum "
            f"{policy.min_rescued_true_positives}"
        )
    if rescued:
        false_positives_per_rescue = added_false_positives / rescued
        if false_positives_per_rescue > policy.max_added_false_positives_per_rescue:
            reasons.append(
                f"added false positives per rescue {false_positives_per_rescue:.3f} exceeds "
                f"{policy.max_added_false_positives_per_rescue:.3f}"
            )
    if union.true_positive < baseline.true_positive:
        reasons.append("union loses baseline true positives")
    if runtime_ratio > policy.max_total_runtime_ratio:
        reasons.append(
            f"runtime ratio {runtime_ratio:.3f} exceeds {policy.max_total_runtime_ratio:.3f}"
        )
    return ToolComparison(
        baseline=baseline.callset,
        candidate=candidate.callset,
        truth_set=baseline.truth_set,
        meets_policy=not reasons,
        rescued_true_positives=rescued,
        added_false_positives=added_false_positives,
        runtime_ratio=round(runtime_ratio, 6),
        reasons=tuple(reasons),
    )


def read_metrics(path: Path) -> BenchmarkMetrics:
    return BenchmarkMetrics.model_validate_json(path.read_text(encoding="utf-8"))


def load_comparison_policy(path: Path) -> ComparisonPolicy:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    return ComparisonPolicy.model_validate(raw)


def write_comparison(comparison: ToolComparison, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(comparison.model_dump_json(indent=2) + "\n", encoding="utf-8")


def metrics_summary(metrics: BenchmarkMetrics) -> str:
    return json.dumps(
        {
            **metrics.model_dump(),
            "precision": round(metrics.precision, 6),
            "recall": round(metrics.recall, 6),
            "f1": round(metrics.f1, 6),
        },
        indent=2,
    )
