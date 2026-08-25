from pathlib import Path


def test_wgs_workflow_uses_samtools_123_supported_flags() -> None:
    workflow = Path("workflow/wgs.smk").read_text(encoding="utf-8")

    for unsupported in ("--sort-by-name", "--name-sort", "--add-mate-score", "--stats"):
        assert unsupported not in workflow
    for supported in ("sort -n", "merge -n", "fixmate -m", "markdup -s"):
        assert supported in workflow
