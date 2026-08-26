"""Typed records crossing the analysis module's single public seam."""

from __future__ import annotations

from enum import StrEnum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

Probability = Annotated[float, Field(ge=0.0, le=1.0)]


class PhaseStatus(StrEnum):
    FAMILY_BACKED_TRANS = "family_backed_trans"
    READ_BACKED_TRANS = "read_backed_trans"
    STATISTICAL_TRANS = "statistical_trans"
    UNRESOLVED = "unresolved"
    INCOMPATIBLE_CIS = "incompatible_cis"


class PhaseMethod(StrEnum):
    FAMILY = "family"
    READ_BACKED = "read_backed"
    STATISTICAL = "statistical"


class VariantKey(BaseModel):
    """A minimal, normalized GRCh38 small-variant identity."""

    model_config = ConfigDict(frozen=True)

    chrom: str
    pos: Annotated[int, Field(ge=1)]
    ref: str
    alt: str

    @field_validator("chrom")
    @classmethod
    def normalize_chrom(cls, value: str) -> str:
        chrom = value.strip()
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        if not chrom:
            raise ValueError("chromosome cannot be empty")
        return chrom.upper() if chrom.upper() in {"X", "Y", "M", "MT"} else chrom

    @field_validator("ref", "alt")
    @classmethod
    def normalize_allele(cls, value: str) -> str:
        allele = value.strip().upper()
        if not allele:
            raise ValueError("allele cannot be empty")
        return allele

    @property
    def challenge_chrom(self) -> str:
        chrom = "M" if self.chrom == "MT" else self.chrom
        return f"chr{chrom}"

    @property
    def label(self) -> str:
        return f"{self.challenge_chrom}:{self.pos}:{self.ref}>{self.alt}"


class SourceReference(BaseModel):
    """One auditable source attached to an evidence assertion."""

    source: str
    release: str
    record_id: str | None = None
    url: str | None = None
    retrieved_at: str | None = None


class GenotypeEvidence(BaseModel):
    """Observed single-sample call evidence."""

    gt: str
    dp: Annotated[int | None, Field(ge=0)] = None
    ad_ref: Annotated[int | None, Field(ge=0)] = None
    ad_alt: Annotated[int | None, Field(ge=0)] = None
    gq: Annotated[float | None, Field(ge=0)] = None
    qual: Annotated[float | None, Field(ge=0)] = None
    filter: str = "PASS"
    phase_set: str | None = None
    phased_gt: str | None = None
    phase_method: PhaseMethod | None = None

    @property
    def allele_balance(self) -> float | None:
        if self.ad_ref is None or self.ad_alt is None:
            return None
        total = self.ad_ref + self.ad_alt
        return self.ad_alt / total if total else None

    @property
    def is_heterozygous(self) -> bool:
        alleles = self.gt.replace("|", "/").split("/")
        return len(alleles) == 2 and set(alleles) == {"0", "1"}

    @property
    def is_homozygous_alt(self) -> bool:
        return self.gt.replace("|", "/") == "1/1"


class TranscriptConsequence(BaseModel):
    """One losslessly retained VEP transcript consequence."""

    transcript: str | None = None
    consequence: str
    biotype: str | None = None
    exon: str | None = None
    intron: str | None = None
    canonical: bool = False
    mane_select: str | None = None
    hgvsc: str | None = None
    hgvsp: str | None = None


class VariantAnnotation(BaseModel):
    """Versioned annotation evidence; missing values remain explicitly unknown."""

    gene: str
    transcript: str | None = None
    consequence: str
    transcripts: tuple[TranscriptConsequence, ...] = ()
    max_population_af: Probability | None = None
    population_ac: Annotated[int | None, Field(ge=0)] = None
    population_an: Annotated[int | None, Field(ge=0)] = None
    clinvar_significance: str | None = None
    clinvar_review_status: str | None = None
    cadd_phred: Annotated[float | None, Field(ge=0)] = None
    revel: Probability | None = None
    spliceai: Probability | None = None
    alphamissense: Probability | None = None
    sift_deleteriousness: Probability | None = None
    polyphen_damagingness: Probability | None = None
    phenotype_gene_score: Probability = 0.0
    disease_mechanism_match: Probability = 0.0
    sources: tuple[SourceReference, ...] = ()


class VariantEvidence(BaseModel):
    key: VariantKey
    genotype: GenotypeEvidence
    annotation: VariantAnnotation


class ScoreContribution(BaseModel):
    term: str
    raw_value: float
    weight: float
    points: float
    rationale: str


class CandidateHypothesis(BaseModel):
    """A ranked causal hypothesis: one variant or an exact same-gene pair."""

    model_config = ConfigDict(frozen=True)

    variants: tuple[VariantEvidence, ...]
    gene: str
    inheritance: Literal["compound_heterozygous", "recessive_homozygous", "single_variant"]
    phase_status: PhaseStatus = PhaseStatus.UNRESOLVED
    score: float
    epcr: Annotated[float, Field(gt=0.0, le=1.0)] = 0.5
    contributions: tuple[ScoreContribution, ...]
    cautions: tuple[str, ...] = ()

    @model_validator(mode="after")
    def check_hypothesis_shape(self) -> CandidateHypothesis:
        expected = 2 if self.inheritance == "compound_heterozygous" else 1
        if len(self.variants) != expected:
            raise ValueError(f"{self.inheritance} requires {expected} variant(s)")
        if any(v.annotation.gene != self.gene for v in self.variants):
            raise ValueError("all variants in a hypothesis must share its gene")
        if len({v.key for v in self.variants}) != len(self.variants):
            raise ValueError("a hypothesis cannot contain the same variant twice")
        return self


class RankedCase(BaseModel):
    proband_id: str
    assembly: Literal["GRCh38"] = "GRCh38"
    policy_version: str
    candidates: tuple[CandidateHypothesis, ...]
    excluded_variant_count: Annotated[int, Field(ge=0)] = 0
    low_support_variant_count: Annotated[int, Field(ge=0)] = 0
    epcr_interpretation: Literal["ordinal"] = "ordinal"
