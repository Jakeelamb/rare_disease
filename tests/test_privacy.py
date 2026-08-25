from pathlib import Path

from mva_hackathon.privacy import audit_public_tree


def test_privacy_audit_detects_token_and_restricted_extension(tmp_path: Path) -> None:
    (tmp_path / "safe.py").write_text("print('safe')\n", encoding="utf-8")
    (tmp_path / "secret.txt").write_text(
        "token=" + "hf_" + "abcdefghijklmnopqrstuvwxyz123456" + "\n", encoding="utf-8"
    )
    (tmp_path / "patient.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

    result = audit_public_tree(
        tmp_path,
        files=[tmp_path / "safe.py", tmp_path / "secret.txt", tmp_path / "patient.vcf"],
    )

    assert not result.ok
    assert {finding.rule for finding in result.findings} == {
        "secret.huggingface_token",
        "restricted.genomic_extension",
    }


def test_privacy_audit_accepts_public_source_tree(tmp_path: Path) -> None:
    source = tmp_path / "module.py"
    source.write_text("VALUE = 'public'\n", encoding="utf-8")

    result = audit_public_tree(tmp_path, files=[source])

    assert result.ok


def test_privacy_audit_rejects_oversized_git_candidate(tmp_path: Path) -> None:
    oversized = tmp_path / "unexpected.bin"
    oversized.write_bytes(b"x" * 17)

    result = audit_public_tree(
        tmp_path,
        files=[oversized],
        max_public_file_bytes=16,
    )

    assert not result.ok
    assert len(result.findings) == 1
    assert result.findings[0].rule == "restricted.oversized_file"
