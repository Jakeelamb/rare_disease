import gzip
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest


def _load_builder() -> ModuleType:
    path = Path("scripts/build_source_reference.py")
    spec = importlib.util.spec_from_file_location("build_source_reference", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_replaces_primary_and_appends_decoys_without_chr_prefix(tmp_path: Path) -> None:
    builder = _load_builder()
    primary = tmp_path / "primary.fa.gz"
    combined = tmp_path / "combined.fa.gz"
    output = tmp_path / "built.fa"
    with gzip.open(primary, "wt", encoding="ascii") as handle:
        handle.write(">chr1 source\nNNGT\n>chr2\nACGT\n")
    with gzip.open(combined, "wt", encoding="ascii") as handle:
        handle.write(">chr1 old\nACGT\n>chr2\nACGT\n>chrDecoy\nTT\n")

    summary = builder.build_reference(primary, combined, output)

    assert summary == {"masked_primary_contigs": 2, "total_contigs": 3}
    assert output.read_text() == ">1 source\nNNGT\n>2\nACGT\n>Decoy\nTT\n"


def test_rejects_primary_order_mismatch(tmp_path: Path) -> None:
    builder = _load_builder()
    primary = tmp_path / "primary.fa"
    combined = tmp_path / "combined.fa"
    output = tmp_path / "built.fa"
    primary.write_text(">chr2\nA\n")
    combined.write_text(">chr1\nA\n")

    with pytest.raises(ValueError, match="primary order mismatch"):
        builder.build_reference(primary, combined, output)
