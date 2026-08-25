"""Union-first comparison of ranked hypotheses with lane provenance."""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import BaseModel, Field

from .models import CandidateHypothesis, RankedCase


class LaneObservation(BaseModel):
    rank: int = Field(ge=1)
    score: float
    epcr: float = Field(gt=0, le=1)


class HypothesisComparison(BaseModel):
    variants: tuple[str, ...]
    gene: str
    lanes: dict[str, LaneObservation]
    review_required: bool


class RankingComparison(BaseModel):
    lane_count: int = Field(ge=2)
    hypothesis_count: int = Field(ge=0)
    consensus_count: int = Field(ge=0)
    hypotheses: tuple[HypothesisComparison, ...]


def _identity(candidate: CandidateHypothesis) -> tuple[str, ...]:
    return tuple(sorted(variant.key.label for variant in candidate.variants))


def compare_ranked_cases(lanes: dict[str, RankedCase], limit: int = 100) -> RankingComparison:
    if len(lanes) < 2:
        raise ValueError("ranking comparison requires at least two named lanes")
    observed: dict[tuple[str, ...], dict[str, tuple[CandidateHypothesis, LaneObservation]]] = {}
    for lane, ranked in sorted(lanes.items()):
        for rank, candidate in enumerate(ranked.candidates[:limit], start=1):
            observed.setdefault(_identity(candidate), {})[lane] = (
                candidate,
                LaneObservation(rank=rank, score=candidate.score, epcr=candidate.epcr),
            )

    rows: list[HypothesisComparison] = []
    for variants, lane_records in observed.items():
        genes = sorted({candidate.gene for candidate, _ in lane_records.values()})
        lane_observations = {
            lane: observation for lane, (_, observation) in sorted(lane_records.items())
        }
        rows.append(
            HypothesisComparison(
                variants=variants,
                gene="|".join(genes),
                lanes=lane_observations,
                review_required=len(lane_records) != len(lanes) or len(genes) != 1,
            )
        )
    rows.sort(
        key=lambda row: (
            min(observation.rank for observation in row.lanes.values()),
            row.gene,
            row.variants,
        )
    )
    consensus = sum(not row.review_required for row in rows)
    return RankingComparison(
        lane_count=len(lanes),
        hypothesis_count=len(rows),
        consensus_count=consensus,
        hypotheses=tuple(rows),
    )


def write_ranking_comparison(comparison: RankingComparison, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".json":
        output.write_text(comparison.model_dump_json(indent=2) + "\n", encoding="utf-8")
        return
    lines = [
        "# Private lane comparison",
        "",
        "> Union-first review queue. Lane agreement is evidence, never a voting threshold.",
        "",
        f"- Lanes: {comparison.lane_count}",
        f"- Union hypotheses: {comparison.hypothesis_count}",
        f"- Present in every lane: {comparison.consensus_count}",
        "",
        "| Gene | Exact variants | Lane ranks | Manual review |",
        "|---|---|---|---|",
    ]
    for row in comparison.hypotheses:
        variants = " + ".join(f"`{item}`" for item in row.variants)
        lane_ranks = ", ".join(
            f"{lane}=#{observation.rank}" for lane, observation in row.lanes.items()
        )
        lines.append(
            f"| {row.gene} | {variants} | {lane_ranks} | "
            f"{'required' if row.review_required else 'concordant'} |"
        )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def comparison_summary(comparison: RankingComparison) -> str:
    return json.dumps(
        {
            "lane_count": comparison.lane_count,
            "hypothesis_count": comparison.hypothesis_count,
            "consensus_count": comparison.consensus_count,
            "review_required_count": sum(item.review_required for item in comparison.hypotheses),
        },
        indent=2,
    )
