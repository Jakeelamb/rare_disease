import hashlib
from pathlib import Path

from mva_hackathon.provenance import build_run_manifest, write_manifest


def test_manifest_hashes_inputs_and_pins_revisions(tmp_path: Path) -> None:
    source = tmp_path / "synthetic.txt"
    source.write_bytes(b"synthetic public fixture\n")

    manifest = build_run_manifest(
        [source],
        {"dataset": "dataset-sha", "space": "space-sha"},
    )
    output = tmp_path / "manifest.json"
    write_manifest(manifest, output)

    assert manifest["challenge_revisions"] == {
        "dataset": "dataset-sha",
        "space": "space-sha",
    }
    assert set(manifest["tools"]) == {
        "bcftools",
        "fastp",
        "multiqc",
        "samtools",
        "snakemake",
    }
    assert manifest["inputs"] == [
        {
            "name": "synthetic.txt",
            "bytes": source.stat().st_size,
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        }
    ]
    assert output.read_text(encoding="utf-8").endswith("\n")
