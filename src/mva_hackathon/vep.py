"""Convert offline Ensembl VEP VCF annotations into typed evidence records."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import yaml
from cyvcf2 import VCF  # type: ignore[import-untyped]

from .models import (
    GenotypeEvidence,
    PhaseMethod,
    SourceReference,
    TranscriptConsequence,
    VariantAnnotation,
    VariantEvidence,
    VariantKey,
)

CSQ_HEADER = re.compile(r"ID=CSQ.*?Format: ([^\"]+)")


def _first_scalar(value: Any) -> Any:
    if value is None:
        return None
    try:
        first = value[0]
        return first[0] if hasattr(first, "__len__") and not isinstance(first, str) else first
    except (IndexError, TypeError):
        return None


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    number = int(value)
    return None if number < 0 else number


def _optional_float(value: Any) -> float | None:
    if value in {None, "", ".", "-"}:
        return None
    number = float(value)
    return None if number < 0 else number


def _format_scalar(record: Any, field: str) -> Any:
    try:
        return _first_scalar(record.format(field))
    except KeyError:
        return None


def _allelic_depths(record: Any) -> tuple[int | None, int | None]:
    try:
        values = record.format("AD")[0]
    except (KeyError, IndexError, TypeError):
        return None, None
    if len(values) < 2:
        return None, None
    return _optional_int(values[0]), _optional_int(values[1])


def _max_frequency(annotation: dict[str, str]) -> float | None:
    values = [
        _optional_float(annotation.get(field))
        for field in ("gnomADe_AF", "gnomADg_AF", "MAX_AF", "AF")
    ]
    present = [value for value in values if value is not None]
    return max(present) if present else None


def load_gene_priors(path: Path) -> tuple[str, dict[str, float]]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    release = str(raw["version"])
    priors = {gene["symbol"]: float(gene["phenotype_gene_score"]) for gene in raw["genes"]}
    return release, priors


def load_mechanism_priors(path: Path) -> tuple[str, dict[str, float]]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    release = str(raw["version"])
    priors = {gene["symbol"]: float(gene["disease_mechanism_match"]) for gene in raw["genes"]}
    return release, priors


def _max_annotation_float(
    annotations: list[dict[str, str]], fields: tuple[str, ...]
) -> float | None:
    values = [
        value
        for annotation in annotations
        for field in fields
        if (value := _optional_float(annotation.get(field))) is not None
    ]
    return max(values) if values else None


def _first_annotation_int(annotations: list[dict[str, str]], fields: tuple[str, ...]) -> int | None:
    for annotation in annotations:
        for field in fields:
            value = annotation.get(field)
            if isinstance(value, str) and value not in {"", ".", "-"}:
                number = int(value)
                if number >= 0:
                    return number
    return None


def _preferred_annotation(annotations: list[dict[str, str]]) -> dict[str, str]:
    return min(
        annotations,
        key=lambda annotation: (
            0 if annotation.get("MANE_SELECT") else 1,
            0 if annotation.get("CANONICAL") == "YES" else 1,
            annotation.get("Feature") or "",
        ),
    )


def vep_to_evidence(
    vcf_path: Path,
    gene_prior_path: Path,
    vep_release: str = "116",
    mechanism_prior_path: Path | None = None,
    phase_method: PhaseMethod | None = None,
) -> list[VariantEvidence]:
    prior_release, phenotype_priors = load_gene_priors(gene_prior_path)
    mechanism_release: str | None = None
    mechanism_priors: dict[str, float] = {}
    if mechanism_prior_path is not None:
        mechanism_release, mechanism_priors = load_mechanism_priors(mechanism_prior_path)
    vcf = VCF(str(vcf_path))
    match = CSQ_HEADER.search(vcf.raw_header)
    if not match:
        raise ValueError("VCF header has no Ensembl VEP CSQ Format declaration")
    fields = match.group(1).split("|")
    results: list[VariantEvidence] = []

    for record in vcf:
        raw_csq = record.INFO.get("CSQ")
        if not raw_csq:
            continue
        genotype = record.genotypes[0]
        left, right, phased = genotype
        gt = f"{left}{'|' if phased else '/'}{right}"
        ad_ref, ad_alt = _allelic_depths(record)
        phase_set = _format_scalar(record, "PID") or _format_scalar(record, "PS")
        phased_gt = _format_scalar(record, "PGT") or (gt if phased else None)

        annotations_by_gene: dict[str, list[dict[str, str]]] = {}
        for item in str(raw_csq).split(","):
            values = item.split("|")
            values.extend([""] * (len(fields) - len(values)))
            annotation = dict(zip(fields, values, strict=True))
            gene = annotation.get("SYMBOL") or annotation.get("Gene")
            consequence = annotation.get("Consequence")
            if not gene or not consequence:
                continue
            annotations_by_gene.setdefault(gene, []).append(annotation)

        for gene, gene_annotations in annotations_by_gene.items():
            preferred = _preferred_annotation(gene_annotations)
            sources: list[SourceReference] = [
                SourceReference(
                    source="Ensembl VEP",
                    release=vep_release,
                    record_id=preferred.get("Feature") or None,
                    url="https://www.ensembl.org/info/docs/tools/vep/",
                ),
                SourceReference(
                    source="Phenotype prior",
                    release=prior_release,
                    record_id=gene,
                ),
            ]
            if mechanism_release is not None:
                sources.append(
                    SourceReference(
                        source="Mechanism prior",
                        release=mechanism_release,
                        record_id=gene,
                    )
                )
            transcript_evidence = tuple(
                TranscriptConsequence(
                    transcript=annotation.get("Feature") or None,
                    consequence=annotation["Consequence"],
                    biotype=annotation.get("BIOTYPE") or None,
                    canonical=annotation.get("CANONICAL") == "YES",
                    mane_select=annotation.get("MANE_SELECT") or None,
                    hgvsc=annotation.get("HGVSc") or None,
                    hgvsp=annotation.get("HGVSp") or None,
                )
                for annotation in gene_annotations
            )
            consequences = sorted(
                {
                    consequence
                    for annotation in gene_annotations
                    for consequence in annotation["Consequence"].split("&")
                    if consequence
                }
            )
            results.append(
                VariantEvidence(
                    key=VariantKey(
                        chrom=str(record.CHROM),
                        pos=int(record.POS),
                        ref=str(record.REF),
                        alt=str(record.ALT[0]),
                    ),
                    genotype=GenotypeEvidence(
                        gt=gt,
                        dp=_optional_int(_format_scalar(record, "DP")),
                        ad_ref=ad_ref,
                        ad_alt=ad_alt,
                        gq=_optional_float(_format_scalar(record, "GQ")),
                        qual=_optional_float(record.QUAL),
                        filter=record.FILTER or "PASS",
                        phase_set=str(phase_set) if phase_set not in {None, ""} else None,
                        phased_gt=str(phased_gt) if phased_gt not in {None, ""} else None,
                        phase_method=phase_method,
                    ),
                    annotation=VariantAnnotation(
                        gene=gene,
                        transcript=preferred.get("Feature") or None,
                        consequence="&".join(consequences),
                        transcripts=transcript_evidence,
                        max_population_af=max(
                            (
                                frequency
                                for annotation in gene_annotations
                                if (frequency := _max_frequency(annotation)) is not None
                            ),
                            default=None,
                        ),
                        population_ac=_first_annotation_int(
                            gene_annotations, ("gnomADe_AC", "gnomADg_AC")
                        ),
                        population_an=_first_annotation_int(
                            gene_annotations, ("gnomADe_AN", "gnomADg_AN")
                        ),
                        clinvar_significance=next(
                            (
                                annotation["CLIN_SIG"]
                                for annotation in gene_annotations
                                if annotation.get("CLIN_SIG")
                            ),
                            None,
                        ),
                        clinvar_review_status=next(
                            (
                                annotation["CLIN_VAR_REVIEW"]
                                for annotation in gene_annotations
                                if annotation.get("CLIN_VAR_REVIEW")
                            ),
                            None,
                        ),
                        cadd_phred=_max_annotation_float(gene_annotations, ("CADD_PHRED",)),
                        revel=_max_annotation_float(gene_annotations, ("REVEL",)),
                        spliceai=_max_annotation_float(
                            gene_annotations,
                            (
                                "SpliceAI_pred_DS_AG",
                                "SpliceAI_pred_DS_AL",
                                "SpliceAI_pred_DS_DG",
                                "SpliceAI_pred_DS_DL",
                            ),
                        ),
                        alphamissense=_max_annotation_float(
                            gene_annotations, ("am_pathogenicity", "AlphaMissense")
                        ),
                        phenotype_gene_score=phenotype_priors.get(gene, 0.0),
                        disease_mechanism_match=mechanism_priors.get(gene, 0.0),
                        sources=tuple(sources),
                    ),
                )
            )
    vcf.close()
    return results


def write_evidence_jsonl(evidence: list[VariantEvidence], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for variant in evidence:
            handle.write(variant.model_dump_json() + "\n")


def read_evidence_jsonl(path: Path) -> list[VariantEvidence]:
    evidence: list[VariantEvidence] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                evidence.append(VariantEvidence.model_validate_json(line))
            except ValueError as error:
                raise ValueError(f"invalid evidence on line {line_number}: {error}") from error
    return evidence


def evidence_summary(evidence: list[VariantEvidence]) -> str:
    genes = {variant.annotation.gene for variant in evidence}
    return json.dumps({"evidence_records": len(evidence), "gene_count": len(genes)}, indent=2)
