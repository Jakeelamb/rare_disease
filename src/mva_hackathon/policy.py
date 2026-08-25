"""Visible scoring-policy schema and loader."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, Field, model_validator


class Thresholds(BaseModel):
    min_depth: int = Field(ge=0)
    min_gq: float = Field(ge=0)
    min_het_allele_balance: float = Field(ge=0, le=1)
    max_het_allele_balance: float = Field(ge=0, le=1)
    max_recessive_af: float = Field(gt=0, le=1)
    max_variants_per_gene: int = Field(ge=2, le=100)
    min_allele_support: float = Field(ge=0, le=1)

    @model_validator(mode="after")
    def validate_balance_range(self) -> Thresholds:
        if self.min_het_allele_balance >= self.max_het_allele_balance:
            raise ValueError("heterozygous allele-balance bounds are reversed")
        return self


class ScoringPolicy(BaseModel):
    version: str
    description: str
    weights: dict[str, float]
    consequence_scores: dict[str, float]
    clinvar_scores: dict[str, float]
    phase_scores: dict[str, float]
    thresholds: Thresholds
    epcr_temperature: float = Field(gt=0)

    @model_validator(mode="after")
    def validate_terms(self) -> ScoringPolicy:
        required = {
            "quality",
            "rarity",
            "consequence",
            "clinical",
            "pathogenicity",
            "phenotype_gene",
            "mechanism",
            "pair",
            "phase",
        }
        missing = required - self.weights.keys()
        if missing:
            raise ValueError(f"missing weights: {sorted(missing)}")
        return self


def load_policy(path: Path) -> ScoringPolicy:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    return ScoringPolicy.model_validate(raw)
