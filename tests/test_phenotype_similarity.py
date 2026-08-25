from pathlib import Path

import pytest

from mva_hackathon.phenotype import PhenotypeManifest, PhenotypeObservation
from mva_hackathon.phenotype_similarity import (
    build_gene_prior_document,
    load_gene_annotations,
    load_ontology,
    score_genes,
)

ONTOLOGY = """format-version: 1.2

[Term]
id: HP:0000001

[Term]
id: HP:0000118
is_a: HP:0000001 ! All

[Term]
id: HP:0001000
alt_id: HP:0099999
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001001
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001002
is_a: HP:0001000 ! Query feature
"""

ANNOTATIONS = """ncbi_gene_id\tgene_symbol\thpo_id\thpo_name
1\tEXACT\tHP:0001000\tQuery feature
2\tDESC\tHP:0001002\tChild feature
3\tSIB\tHP:0001001\tSibling feature
"""


def _write_sources(tmp_path: Path) -> tuple[Path, Path]:
    ontology_path = tmp_path / "hp.obo"
    annotation_path = tmp_path / "genes_to_phenotype.txt"
    ontology_path.write_text(ONTOLOGY, encoding="utf-8")
    annotation_path.write_text(ANNOTATIONS, encoding="utf-8")
    return ontology_path, annotation_path


def test_query_directed_resnik_rewards_exact_and_descendant_terms(tmp_path: Path) -> None:
    ontology_path, annotation_path = _write_sources(tmp_path)
    ontology = load_ontology(ontology_path)
    annotations = load_gene_annotations(annotation_path, ontology)

    scores = score_genes({"HP:0099999"}, annotations, ontology)

    assert scores["EXACT"] == 1.0
    assert scores["DESC"] == 1.0
    assert scores["SIB"] == 0.0


def test_prior_document_excludes_family_terms_and_mechanism_knowledge(tmp_path: Path) -> None:
    ontology_path, annotation_path = _write_sources(tmp_path)
    phenotype = PhenotypeManifest(
        source_name="private.docx",
        source_sha256="abc123",
        curation_status="reviewed",
        reviewer="synthetic-test",
        reviewed_at="2026-08-25T00:00:00+00:00",
        observations=(
            PhenotypeObservation(
                feature="query",
                hpo_term="query",
                hpo_id="HP:0001000",
                notes="",
                subject="proband",
            ),
            PhenotypeObservation(
                feature="relative",
                hpo_term="sibling",
                hpo_id="HP:0001001",
                notes="",
                subject="family",
            ),
        ),
    )

    document = build_gene_prior_document(
        phenotype,
        ontology_path,
        annotation_path,
        "test-hpo-v1",
    )
    genes = {gene["symbol"]: gene for gene in document["genes"]}

    assert document["version"] == "hpo-test-hpo-v1"
    assert document["sources"]["phenotype_source_sha256"] == "abc123"
    assert genes["EXACT"]["phenotype_gene_score"] == 1.0
    assert genes["SIB"]["phenotype_gene_score"] == 0.0
    assert all("disease_mechanism_match" not in gene for gene in document["genes"])


def test_gene_scoring_rejects_unreviewed_automatic_extraction(tmp_path: Path) -> None:
    ontology_path, annotation_path = _write_sources(tmp_path)
    phenotype = PhenotypeManifest(
        source_name="private.docx",
        source_sha256="abc123",
        observations=(
            PhenotypeObservation(
                feature="query",
                hpo_term="query",
                hpo_id="HP:0001000",
                notes="",
                subject="unresolved",
            ),
        ),
    )

    with pytest.raises(ValueError, match="manual phenotype curation"):
        build_gene_prior_document(phenotype, ontology_path, annotation_path, "test-hpo-v1")
