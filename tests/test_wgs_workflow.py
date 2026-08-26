from pathlib import Path


def test_wgs_workflow_uses_samtools_123_supported_flags() -> None:
    workflow = Path("workflow/wgs.smk").read_text(encoding="utf-8")

    for unsupported in ("--sort-by-name", "--name-sort", "--add-mate-score", "--stats"):
        assert unsupported not in workflow
    for supported in ("sort -n", "merge -n", "fixmate -m", "markdup -s"):
        assert supported in workflow


def test_all_vep_lanes_emit_reference_backed_hgvs_and_native_missense_scores() -> None:
    workflow = Path("workflow/Snakefile").read_text(encoding="utf-8")

    assert workflow.count("--fasta {input.reference:q}") == 3
    assert workflow.count("--hgvs --sift b --polyphen b") == 3
    assert workflow.count("--plugin AlphaMissense,file={input.alphamissense:q}") == 3
    assert workflow.count("--phase-method read_backed") == 3


def test_wgs_comparison_ends_in_local_manual_adjudication_not_submission() -> None:
    workflow = Path("workflow/Snakefile").read_text(encoding="utf-8")
    wgs_workflow = Path("workflow/wgs.smk").read_text(encoding="utf-8")

    assert "query-candidate-clinvar" in workflow
    assert "adjudicate-leading" in workflow
    assert "agnostic_top30_reads.review.md" in wgs_workflow
    assert "leading_candidate.review.md" in wgs_workflow
    assert "create-submission" not in workflow
    assert "upload" not in wgs_workflow.lower()


def test_direct_read_phase_threshold_is_explicit_in_workflow() -> None:
    workflow = Path("workflow/wgs.smk").read_text(encoding="utf-8")

    assert workflow.count("--min-phase-support-fragments") == 2


def test_wgs_phasing_scatters_primary_contigs_and_preserves_other_calls() -> None:
    workflow = Path("workflow/wgs.smk").read_text(encoding="utf-8")

    for rule in (
        "rule split_deepvariant_phase_contig:",
        "rule phase_deepvariant_contig:",
        "rule retain_unphased_nonprimary_calls:",
        "rule phase_deepvariant:",
    ):
        assert rule in workflow
    assert "bcftools" in workflow
    assert " concat --threads " in workflow
    assert "mem_mb=12000" in workflow
