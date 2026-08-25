import pytest

from mva_hackathon.phenotype import (
    PhenotypeManifest,
    PhenotypeObservation,
    curate_phenotype,
    normalize_hpo_id,
    validate_phenotype_curation,
)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("HP:0002859", "HP:0002859"),
        ("HP:000285\n9", "HP:0002859"),
        (" hp:0000121 ", "HP:0000121"),
    ],
)
def test_normalize_hpo_id_handles_docx_line_wrapping(raw: str, expected: str) -> None:
    assert normalize_hpo_id(raw) == expected


def test_normalize_hpo_id_rejects_malformed_values() -> None:
    with pytest.raises(ValueError, match="no valid HPO"):
        normalize_hpo_id("HP:123")


def test_manual_curation_requires_an_explicit_decision_for_every_row() -> None:
    raw = PhenotypeManifest(
        source_name="synthetic.docx",
        source_sha256="abc",
        observations=(
            PhenotypeObservation(
                feature="feature one",
                hpo_term="one",
                hpo_id="HP:0000001",
                notes="",
                subject="unresolved",
            ),
            PhenotypeObservation(
                feature="feature two",
                hpo_term="two",
                hpo_id="HP:0000002",
                notes="",
                subject="unresolved",
            ),
        ),
    )

    with pytest.raises(ValueError, match="exactly one decision for every phenotype row"):
        curate_phenotype(raw, {1: "proband_present"}, reviewer="Jake")

    curated = curate_phenotype(
        raw,
        {1: "proband_present", 2: "family_absent"},
        reviewer="Jake",
    )

    validate_phenotype_curation(curated)
    assert curated.curation_status == "reviewed"
    assert curated.observations[0].subject == "proband"
    assert curated.observations[0].present
    assert curated.observations[1].subject == "family"
    assert not curated.observations[1].present
