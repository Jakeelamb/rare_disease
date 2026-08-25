import pytest

from mva_hackathon.benchmark import (
    BenchmarkMetrics,
    ComparisonPolicy,
    compare_tool_metrics,
)


def metrics(name: str, *, tp: int, fp: int, fn: int, runtime: float = 100.0) -> BenchmarkMetrics:
    return BenchmarkMetrics(
        callset=name,
        truth_set="HG002-GRCh38",
        engine="vcfeval-3.13",
        true_positive=tp,
        false_positive=fp,
        false_negative=fn,
        runtime_seconds=runtime,
        peak_rss_mb=1024,
    )


def test_reports_bounded_incremental_rescue() -> None:
    comparison = compare_tool_metrics(
        baseline=metrics("supplied", tp=980, fp=10, fn=20),
        candidate=metrics("deepvariant", tp=985, fp=12, fn=15),
        union=metrics("union", tp=990, fp=14, fn=10, runtime=210),
        policy=ComparisonPolicy(
            min_rescued_true_positives=2,
            max_added_false_positives_per_rescue=1.0,
            max_total_runtime_ratio=3.0,
        ),
    )

    assert comparison.meets_policy
    assert comparison.rescued_true_positives == 10
    assert comparison.added_false_positives == 4


def test_reports_when_an_extra_tool_adds_no_value() -> None:
    comparison = compare_tool_metrics(
        baseline=metrics("supplied", tp=980, fp=10, fn=20),
        candidate=metrics("redundant", tp=979, fp=30, fn=21),
        union=metrics("union", tp=980, fp=28, fn=20, runtime=250),
        policy=ComparisonPolicy(),
    )

    assert not comparison.meets_policy
    assert any("no incremental true-positive rescue" in reason for reason in comparison.reasons)


def test_metrics_expose_derived_values_without_claiming_calibration() -> None:
    observed = metrics("caller", tp=80, fp=20, fn=20)

    assert observed.precision == 0.8
    assert observed.recall == 0.8
    assert observed.f1 == pytest.approx(0.8)
