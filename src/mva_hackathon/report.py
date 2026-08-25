"""Deterministic private evidence-card rendering for human review."""

from __future__ import annotations

from pathlib import Path

from .models import RankedCase


def render_ranked_case(ranked: RankedCase, limit: int = 10) -> str:
    lines = [
        "# Private candidate review",
        "",
        "> Research prioritization only. Not a diagnosis or treatment recommendation.",
        "",
        f"- Proband key: `{ranked.proband_id}`",
        f"- Assembly: `{ranked.assembly}`",
        f"- Policy: `{ranked.policy_version}`",
        f"- Eligible hypotheses: {len(ranked.candidates)}",
        f"- Excluded input variants: {ranked.excluded_variant_count}",
        f"- Low-support variants: {ranked.low_support_variant_count}",
        "- EPCR interpretation: ordinal ordering heuristic; not a calibrated probability",
        "",
    ]
    for rank, candidate in enumerate(ranked.candidates[:limit], start=1):
        alleles = " + ".join(f"`{variant.key.label}`" for variant in candidate.variants)
        lines.extend(
            [
                f"## {rank}. {candidate.gene}",
                "",
                f"- Hypothesis: `{candidate.inheritance}`",
                f"- Alleles: {alleles}",
                f"- Score / ordinal EPCR: `{candidate.score:.6f}` / `{candidate.epcr:.9f}`",
                f"- Phase: `{candidate.phase_status.value}`",
                f"- Cautions: {'; '.join(candidate.cautions) or 'none recorded'}",
                "",
                "| Term | Raw | Weight | Points | Rationale |",
                "|---|---:|---:|---:|---|",
            ]
        )
        for item in candidate.contributions:
            rationale = item.rationale.replace("|", "\\|").replace("\n", " ")
            lines.append(
                f"| {item.term} | {item.raw_value:.6f} | {item.weight:.3f} | "
                f"{item.points:.6f} | {rationale} |"
            )
        lines.append("")
        lines.append("Evidence sources:")
        lines.append("")
        sources = {
            (source.source, source.release, source.record_id or "")
            for variant in candidate.variants
            for source in variant.annotation.sources
        }
        if sources:
            lines.extend(
                f"- {source} release `{release}` record `{record_id or 'n/a'}`"
                for source, release, record_id in sorted(sources)
            )
        else:
            lines.append("- No source records attached.")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_ranked_case_report(ranked: RankedCase, output: Path, limit: int = 10) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_ranked_case(ranked, limit=limit), encoding="utf-8")
