from pathlib import Path

import yaml

from mva_hackathon.models import PhaseMethod
from mva_hackathon.vep import read_evidence_jsonl, vep_to_evidence, write_evidence_jsonl


def _phenotype_priors(tmp_path: Path) -> Path:
    path = tmp_path / "phenotype-priors.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "version": "hpo-test-v1",
                "genes": [{"symbol": "BUB1B", "phenotype_gene_score": 1.0}],
            }
        ),
        encoding="utf-8",
    )
    return path


def test_vep_parser_preserves_all_transcripts_predictions_and_phenotype_priors(
    tmp_path: Path,
) -> None:
    evidence = vep_to_evidence(
        Path("tests/fixtures/vep.synthetic.vcf.txt"), _phenotype_priors(tmp_path)
    )

    assert len(evidence) == 2
    assert {item.annotation.gene for item in evidence} == {"BUB1B"}
    assert evidence[0].annotation.max_population_af == 0.00002
    assert evidence[0].annotation.phenotype_gene_score == 1.0
    assert evidence[0].annotation.disease_mechanism_match == 0.0
    assert len(evidence[0].annotation.transcripts) == 2
    assert all(item.exon is None for item in evidence[0].annotation.transcripts)
    assert {item.transcript for item in evidence[0].annotation.transcripts} == {
        "ENST1",
        "ENST2",
    }
    assert evidence[0].annotation.revel == 0.82
    assert evidence[0].annotation.spliceai == 0.71
    assert evidence[0].annotation.sift_deleteriousness == 0.99
    assert evidence[0].annotation.polyphen_damagingness == 0.98
    assert evidence[0].annotation.alphamissense == 0.77
    assert any(source.source == "AlphaMissense" for source in evidence[0].annotation.sources)
    assert evidence[0].genotype.allele_balance == 21 / 41
    assert evidence[0].genotype.phase_set == "block1"
    assert evidence[0].genotype.phased_gt == "0|1"
    assert evidence[0].genotype.phase_method is None

    output = tmp_path / "evidence.jsonl"
    write_evidence_jsonl(evidence, output)
    assert read_evidence_jsonl(output) == evidence


def test_physical_phase_metadata_requires_an_explicit_evidence_method(tmp_path: Path) -> None:
    evidence = vep_to_evidence(
        Path("tests/fixtures/vep.synthetic.vcf.txt"),
        _phenotype_priors(tmp_path),
        phase_method=PhaseMethod.READ_BACKED,
    )

    assert evidence[0].genotype.phase_set == "block1"
    assert evidence[0].genotype.phased_gt == "0|1"
    assert evidence[0].genotype.phase_method == PhaseMethod.READ_BACKED
    assert evidence[1].genotype.phase_set == "block1"
    assert evidence[1].genotype.phased_gt == "1|0"
    assert evidence[1].genotype.phase_method == PhaseMethod.READ_BACKED


def test_mechanism_prior_is_explicit_opt_in_lane(tmp_path: Path) -> None:
    mechanism = tmp_path / "mechanism.yaml"
    mechanism.write_text(
        yaml.safe_dump(
            {
                "version": "mechanism-test-v1",
                "genes": [{"symbol": "BUB1B", "disease_mechanism_match": 1.0}],
            }
        ),
        encoding="utf-8",
    )

    agnostic = vep_to_evidence(
        Path("tests/fixtures/vep.synthetic.vcf.txt"), _phenotype_priors(tmp_path)
    )
    informed = vep_to_evidence(
        Path("tests/fixtures/vep.synthetic.vcf.txt"),
        _phenotype_priors(tmp_path),
        mechanism_prior_path=mechanism,
    )

    assert all(item.annotation.disease_mechanism_match == 0.0 for item in agnostic)
    assert all(item.annotation.disease_mechanism_match == 1.0 for item in informed)
