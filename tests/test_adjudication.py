import csv
import gzip
from pathlib import Path

import pysam

from mva_hackathon.adjudication import (
    ClinVarCandidateQuery,
    VcfCallObservation,
    build_leading_adjudication,
    inspect_exact_vcf_calls,
    query_clinvar_archive,
    render_leading_adjudication,
)
from mva_hackathon.models import (
    CandidateHypothesis,
    GenotypeEvidence,
    RankedCase,
    TranscriptConsequence,
    VariantAnnotation,
    VariantEvidence,
    VariantKey,
)
from mva_hackathon.read_evidence import (
    AlleleReadEvidence,
    CandidateReadAudit,
    PairReadEvidence,
)


def _variant(position: int, ref: str, alt: str, consequence: str) -> VariantEvidence:
    return VariantEvidence(
        key=VariantKey(chrom="1", pos=position, ref=ref, alt=alt),
        genotype=GenotypeEvidence(
            gt="0/1",
            dp=40,
            ad_ref=20,
            ad_alt=20,
            gq=99,
        ),
        annotation=VariantAnnotation(
            gene="SYNTHETIC",
            consequence=consequence,
            transcripts=(
                TranscriptConsequence(
                    transcript="ENST_TEST",
                    consequence=consequence,
                    mane_select="NM_TEST.1",
                    hgvsc=f"ENST_TEST:c.{position}A>G",
                    hgvsp="ENSP_TEST:p.Test",
                ),
            ),
            max_population_af=0.00001,
            alphamissense=0.91 if consequence == "missense_variant" else None,
            sift_deleteriousness=0.98 if consequence == "missense_variant" else None,
            polyphen_damagingness=0.99 if consequence == "missense_variant" else None,
            phenotype_gene_score=0.9,
            disease_mechanism_match=1.0,
        ),
    )


def _ranked_pair() -> RankedCase:
    candidate = CandidateHypothesis(
        variants=(
            _variant(100, "A", "G", "stop_gained"),
            _variant(110, "C", "T", "missense_variant"),
        ),
        gene="SYNTHETIC",
        inheritance="compound_heterozygous",
        score=50,
        epcr=0.8,
        contributions=(),
        cautions=("trans phase is unresolved",),
    )
    return RankedCase(proband_id="PROBAND01", policy_version="test", candidates=(candidate,))


def _write_vcf(tmp_path: Path) -> Path:
    raw = tmp_path / "calls.vcf"
    raw.write_text(
        """##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND01
1\t100\t.\tA\tG\t60\tPASS\t.\tGT:DP:AD:GQ:PS\t0|1:42:21,21:99:100
1\t110\t.\tC\tT\t70\tPASS\t.\tGT:DP:AD:GQ:PS\t1|0:38:18,20:98:100
1\t120\t.\tA\tC,G\t80\tPASS\t.\tGT:DP:AD:GQ:PS\t0/2:20:10,0,10:97:.
""",
        encoding="utf-8",
    )
    compressed = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(raw), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def _write_clinvar(tmp_path: Path) -> Path:
    archive = tmp_path / "variant_summary.txt.gz"
    fields = [
        "#AlleleID",
        "Name",
        "GeneSymbol",
        "ClinicalSignificance",
        "LastEvaluated",
        "RS# (dbSNP)",
        "RCVaccession",
        "PhenotypeList",
        "Assembly",
        "Chromosome",
        "ReviewStatus",
        "NumberSubmitters",
        "VariationID",
        "PositionVCF",
        "ReferenceAlleleVCF",
        "AlternateAlleleVCF",
    ]
    rows = [
        {
            "#AlleleID": "10",
            "Name": "exact nucleotide",
            "GeneSymbol": "SYNTHETIC",
            "ClinicalSignificance": "Pathogenic",
            "LastEvaluated": "Jan 01, 2026",
            "RS# (dbSNP)": "123",
            "RCVaccession": "RCV0001",
            "PhenotypeList": "Synthetic condition",
            "Assembly": "GRCh38",
            "Chromosome": "1",
            "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
            "NumberSubmitters": "2",
            "VariationID": "20",
            "PositionVCF": "100",
            "ReferenceAlleleVCF": "A",
            "AlternateAlleleVCF": "G",
        },
        {
            "#AlleleID": "11",
            "Name": "different nucleotide same protein consequence",
            "GeneSymbol": "SYNTHETIC",
            "ClinicalSignificance": "Uncertain significance",
            "LastEvaluated": "Jan 02, 2026",
            "RS# (dbSNP)": "-",
            "RCVaccession": "RCV0002",
            "PhenotypeList": "Synthetic condition",
            "Assembly": "GRCh38",
            "Chromosome": "1",
            "ReviewStatus": "criteria provided, single submitter",
            "NumberSubmitters": "1",
            "VariationID": "21",
            "PositionVCF": "100",
            "ReferenceAlleleVCF": "A",
            "AlternateAlleleVCF": "T",
        },
    ]
    with gzip.open(archive, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return archive


def test_exact_vcf_recall_preserves_quality_and_read_backed_trans(tmp_path: Path) -> None:
    keys = tuple(variant.key for variant in _ranked_pair().candidates[0].variants)

    calls = inspect_exact_vcf_calls(_write_vcf(tmp_path), keys, caller="independent")

    assert [call.present for call in calls] == [True, True]
    assert [(call.gt, call.dp, call.ad_ref, call.ad_alt) for call in calls] == [
        ("0|1", 42, 21, 21),
        ("1|0", 38, 18, 20),
    ]
    assert [call.alt_haplotype for call in calls] == [1, 0]


def test_exact_vcf_alt_record_is_not_recalled_when_genotype_selects_another_alt(
    tmp_path: Path,
) -> None:
    candidate = VariantKey(chrom="1", pos=120, ref="A", alt="C")

    (call,) = inspect_exact_vcf_calls(_write_vcf(tmp_path), (candidate,), caller="independent")

    assert call.sample == "PROBAND01"
    assert call.gt == "0/2"
    assert call.ad_alt == 0
    assert not call.present


def test_clinvar_query_requires_exact_nucleotide_identity(tmp_path: Path) -> None:
    keys = tuple(variant.key for variant in _ranked_pair().candidates[0].variants)

    query = query_clinvar_archive(keys, _write_clinvar(tmp_path), release="test-2026")

    assert query.queried_variants == keys
    assert len(query.matches) == 1
    assert query.matches[0].variant == keys[0]
    assert query.matches[0].clinical_significance == "Pathogenic"
    assert query.matches[0].number_submitters == 2


def test_adjudication_integrates_ranks_calls_reads_and_exact_clinvar(tmp_path: Path) -> None:
    focused = _ranked_pair()
    keys = tuple(variant.key for variant in focused.candidates[0].variants)
    calls = inspect_exact_vcf_calls(_write_vcf(tmp_path), keys, caller="DeepVariant+WhatsHap")
    clinvar = query_clinvar_archive(keys, _write_clinvar(tmp_path), release="test-2026")
    reads = CandidateReadAudit(
        bam="private.bam",
        min_mapping_quality=20,
        min_base_quality=20,
        alleles=tuple(
            AlleleReadEvidence(
                variant=key,
                usable_depth=40,
                ref_reads=20,
                alt_reads=20,
                other_reads=0,
                filtered_reads=0,
                alt_forward_reads=10,
                alt_reverse_reads=10,
                mean_alt_mapping_quality=60,
                mean_alt_base_quality=39,
            )
            for key in keys
        ),
        pairs=(
            PairReadEvidence(
                left=keys[0],
                right=keys[1],
                shared_fragments=0,
                alt_alt_fragments=0,
                alt_ref_fragments=0,
                ref_alt_fragments=0,
                ref_ref_fragments=0,
                other_fragments=0,
            ),
        ),
    )

    adjudication = build_leading_adjudication(
        focused,
        {"focused": focused, "agnostic": focused, "recall": focused},
        calls,
        reads,
        ClinVarCandidateQuery.model_validate(clinvar),
    )
    report = render_leading_adjudication(adjudication)

    assert adjudication.recall_status == "both_exact_alleles_recalled"
    assert adjudication.recall_pair_phase == "read_backed_trans"
    assert adjudication.lane_ranks == {"agnostic": 1, "focused": 1, "recall": 1}
    assert adjudication.pair_read_evidence is not None
    assert adjudication.pair_read_evidence.interpretation == "no_spanning_fragments"
    assert "no automatic ACMG/AMP classification" in report
    assert "different nucleotide" not in report
    assert "Caller concordance uses the same reads and reference" in report


def test_adjudication_does_not_infer_phase_without_shared_block() -> None:
    focused = _ranked_pair()
    keys = tuple(variant.key for variant in focused.candidates[0].variants)
    calls = tuple(
        VcfCallObservation(caller="recall", variant=key, present=True, gt="0/1") for key in keys
    )
    empty_audit = CandidateReadAudit(
        bam="private.bam",
        min_mapping_quality=20,
        min_base_quality=20,
        alleles=(),
        pairs=(),
    )
    empty_clinvar = ClinVarCandidateQuery(
        release="test",
        archive="archive.gz",
        queried_variants=keys,
        matches=(),
    )

    adjudication = build_leading_adjudication(
        focused,
        {"focused": focused},
        calls,
        empty_audit,
        empty_clinvar,
    )

    assert adjudication.recall_pair_phase == "no_shared_phase_block"
