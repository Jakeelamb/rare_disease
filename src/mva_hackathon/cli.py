"""Command-line interface for auditable workflow entry points."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from .adjudication import (
    ClinVarCandidateQuery,
    build_leading_adjudication,
    inspect_exact_vcf_calls,
    leading_variant_keys,
    query_clinvar_archive,
    write_clinvar_query,
    write_leading_adjudication,
)
from .benchmark import (
    compare_tool_metrics,
    load_comparison_policy,
    read_metrics,
    write_comparison,
)
from .comparison import (
    compare_ranked_cases,
    comparison_summary,
    write_ranking_comparison,
)
from .models import PhaseMethod, RankedCase
from .phenotype import (
    PhenotypeManifest,
    curate_phenotype,
    extract_phenotype,
    public_summary,
    render_phenotype_review,
    validate_phenotype_curation,
    write_phenotype,
)
from .phenotype_similarity import build_gene_prior_document, write_gene_prior_document
from .policy import load_policy
from .privacy import audit_public_tree
from .provenance import build_run_manifest, write_manifest
from .ranking import rank_case
from .read_evidence import CandidateReadAudit, inspect_ranked_bam, write_candidate_read_audit
from .report import write_ranked_case_report
from .submission import validate_submission
from .submission_package import create_submission_package, sha256_file, validate_methods_report
from .vcf import inspect_vcf
from .vep import evidence_summary, read_evidence_jsonl, vep_to_evidence, write_evidence_jsonl

app = typer.Typer(no_args_is_help=True, pretty_exceptions_show_locals=False)
console = Console()


@app.command("privacy-audit")
def privacy_audit_command(
    root: Annotated[Path, typer.Option("--root", file_okay=False, dir_okay=True)] = Path("."),
) -> None:
    """Fail if public Git candidates contain patient-capable files or credentials."""
    result = audit_public_tree(root)
    for finding in result.findings:
        console.print(f"[red]ERROR[/red] {finding.rule} {finding.path}: {finding.detail}")
    if not result.ok:
        raise typer.Exit(code=2)
    console.print("[green]PASS[/green] public-tree privacy audit")


@app.command("compare-rankings")
def compare_rankings_command(
    rankings: Annotated[list[str], typer.Argument(help="Repeated LABEL=RANKED_JSON inputs")],
    output: Annotated[Path, typer.Option("--output", "-o")],
    limit: Annotated[int, typer.Option("--limit", min=1)] = 100,
) -> None:
    """Build a private union-first discrepancy queue across analysis lanes."""
    lanes: dict[str, RankedCase] = {}
    for item in rankings:
        label, separator, raw_path = item.partition("=")
        path = Path(raw_path)
        if not separator or not label or not path.is_file():
            raise typer.BadParameter(f"ranking must be LABEL=existing.json: {item}")
        if label in lanes:
            raise typer.BadParameter(f"duplicate ranking label: {label}")
        lanes[label] = RankedCase.model_validate_json(path.read_text(encoding="utf-8"))
    comparison = compare_ranked_cases(lanes, limit=limit)
    write_ranking_comparison(comparison, output)
    console.print(comparison_summary(comparison))


@app.command("compare-benchmarks")
def compare_benchmarks_command(
    baseline: Annotated[Path, typer.Option("--baseline", exists=True, dir_okay=False)],
    candidate: Annotated[Path, typer.Option("--candidate", exists=True, dir_okay=False)],
    union: Annotated[Path, typer.Option("--union", exists=True, dir_okay=False)],
    policy: Annotated[Path, typer.Option("--policy", exists=True, dir_okay=False)],
    output: Annotated[Path, typer.Option("--output", "-o")],
) -> None:
    """Compare a candidate tool with a baseline on public truth data."""
    comparison = compare_tool_metrics(
        baseline=read_metrics(baseline),
        candidate=read_metrics(candidate),
        union=read_metrics(union),
        policy=load_comparison_policy(policy),
    )
    write_comparison(comparison, output)
    console.print(comparison.model_dump_json(indent=2))


@app.command("manifest")
def manifest_command(
    inputs: Annotated[list[Path], typer.Argument(exists=True, dir_okay=False, readable=True)],
    output: Annotated[Path, typer.Option("--output", "-o")],
    dataset_revision: Annotated[str, typer.Option("--dataset-revision")],
    space_revision: Annotated[str, typer.Option("--space-revision")],
) -> None:
    """Hash immutable inputs and write a machine-readable local provenance manifest."""
    manifest = build_run_manifest(
        inputs,
        {
            "dataset": dataset_revision,
            "space": space_revision,
        },
    )
    write_manifest(manifest, output)
    console.print(
        json.dumps(
            {
                "hashed_inputs": len(inputs),
                "total_bytes": sum(item["bytes"] for item in manifest["inputs"]),
                "manifest": str(output),
            },
            indent=2,
        )
    )


@app.command("extract-phenotype")
def extract_phenotype_command(
    document: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    output: Annotated[Path, typer.Option("--output", "-o")],
) -> None:
    """Extract the restricted DOCX locally; print only a non-identifying summary."""
    manifest = extract_phenotype(document)
    write_phenotype(manifest, output)
    console.print(public_summary(manifest))


@app.command("render-phenotype-review")
def render_phenotype_review_command(
    phenotype_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    output: Annotated[Path, typer.Option("--output", "-o")],
) -> None:
    """Render the private row-indexed curation worksheet."""
    manifest = PhenotypeManifest.model_validate_json(phenotype_json.read_text(encoding="utf-8"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_phenotype_review(manifest), encoding="utf-8")
    console.print(json.dumps({"review": str(output), "rows": len(manifest.observations)}))


@app.command("curate-phenotype")
def curate_phenotype_command(
    phenotype_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    decisions: Annotated[
        list[str], typer.Option("--decision", help="ROW=DECISION; repeat for every row")
    ],
    reviewer: Annotated[str, typer.Option("--reviewer")],
    output: Annotated[Path, typer.Option("--output", "-o")],
) -> None:
    """Create a reviewed phenotype manifest from explicit row decisions."""
    parsed: dict[int, str] = {}
    for item in decisions:
        raw_row, separator, decision = item.partition("=")
        if not separator or not raw_row.isdigit() or int(raw_row) in parsed:
            raise typer.BadParameter(f"decision must be a unique ROW=DECISION: {item}")
        parsed[int(raw_row)] = decision
    raw = PhenotypeManifest.model_validate_json(phenotype_json.read_text(encoding="utf-8"))
    curated = curate_phenotype(raw, parsed, reviewer=reviewer)
    write_phenotype(curated, output)
    console.print(
        json.dumps(
            {"curated_manifest": str(output), "reviewed_rows": len(curated.observations)},
            indent=2,
        )
    )


@app.command("validate-phenotype-curation")
def validate_phenotype_curation_command(
    phenotype_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
) -> None:
    """Fail unless every extracted phenotype row has an explicit manual decision."""
    manifest = PhenotypeManifest.model_validate_json(phenotype_json.read_text(encoding="utf-8"))
    validate_phenotype_curation(manifest)
    console.print("[green]PASS[/green] phenotype curation is complete")


@app.command("inspect-vcf")
def inspect_vcf_command(
    vcf: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
) -> None:
    """Print aggregate VCF structure/QC only; never print patient variant rows."""
    console.print_json(json.dumps(inspect_vcf(vcf)))


@app.command("phenotype-priors")
def phenotype_priors_command(
    phenotype_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    ontology: Annotated[Path, typer.Option(exists=True, dir_okay=False, readable=True)],
    annotations: Annotated[Path, typer.Option(exists=True, dir_okay=False, readable=True)],
    hpo_release: Annotated[str, typer.Option("--hpo-release")],
    output: Annotated[Path, typer.Option("--output", "-o")],
) -> None:
    """Build local, source-versioned HPO gene scores without hosted APIs."""
    phenotype = PhenotypeManifest.model_validate_json(phenotype_json.read_text(encoding="utf-8"))
    document = build_gene_prior_document(
        phenotype,
        ontology,
        annotations,
        hpo_release,
    )
    write_gene_prior_document(document, output)
    console.print(
        json.dumps(
            {
                "gene_count": len(document["genes"]),
                "hpo_release": hpo_release,
                "output": str(output),
            },
            indent=2,
        )
    )


@app.command("vep-to-evidence")
def vep_to_evidence_command(
    vcf: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    gene_priors: Annotated[
        Path, typer.Option("--gene-priors", exists=True, dir_okay=False, readable=True)
    ],
    output: Annotated[Path, typer.Option("--output", "-o")],
    mechanism_priors: Annotated[
        Path | None,
        typer.Option("--mechanism-priors", exists=True, dir_okay=False, readable=True),
    ] = None,
    phase_method: Annotated[PhaseMethod | None, typer.Option("--phase-method")] = None,
) -> None:
    """Convert offline VEP annotations to local typed evidence JSONL."""
    evidence = vep_to_evidence(
        vcf,
        gene_priors,
        mechanism_prior_path=mechanism_priors,
        phase_method=phase_method,
    )
    write_evidence_jsonl(evidence, output)
    console.print(evidence_summary(evidence))


@app.command("rank")
def rank_command(
    evidence_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    policy_path: Annotated[
        Path, typer.Option("--policy", exists=True, dir_okay=False, readable=True)
    ],
    ranked_json: Annotated[Path, typer.Option("--ranked-json")],
) -> None:
    """Rank typed evidence for research review."""
    evidence = read_evidence_jsonl(evidence_path)
    ranked = rank_case(evidence, load_policy(policy_path))
    ranked_json.parent.mkdir(parents=True, exist_ok=True)
    ranked_json.write_text(ranked.model_dump_json(indent=2) + "\n", encoding="utf-8")
    console.print(
        json.dumps(
            {
                "input_variants": len(evidence),
                "ranked_hypotheses": len(ranked.candidates),
                "policy_version": ranked.policy_version,
            },
            indent=2,
        )
    )


@app.command("validate-submission")
def validate_submission_command(
    csv_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
) -> None:
    """Validate the exact public Track 1 CSV contract before using an upload."""
    result = validate_submission(csv_path)
    for warning in result.warnings:
        console.print(f"[yellow]WARNING[/yellow] {warning}")
    for error in result.errors:
        console.print(f"[red]ERROR[/red] {error}")
    if not result.ok:
        raise typer.Exit(code=2)
    console.print("[green]PASS[/green] submission contract is valid")


@app.command("submission-preflight")
def submission_preflight_command(
    csv_path: Annotated[
        Path, typer.Option("--csv", exists=True, dir_okay=False, readable=True)
    ] = Path("results/private/submission/final/jakeelamb_jl_wgs_v1.csv"),
    report_path: Annotated[
        Path, typer.Option("--report", exists=True, dir_okay=False, readable=True)
    ] = Path("results/private/submission/final/jakeelamb_track1_report.md"),
) -> None:
    """Validate reviewed CSV/report artifacts without creating or uploading them."""
    csv_result = validate_submission(csv_path)
    report_result = validate_methods_report(report_path)
    payload = {
        "csv": {
            "path": str(csv_path),
            "sha256": sha256_file(csv_path),
            "status": "PASS" if csv_result.ok else "FAIL",
            "errors": list(csv_result.errors),
            "warnings": list(csv_result.warnings),
        },
        "report": {
            "path": str(report_path),
            "sha256": sha256_file(report_path),
            "status": "PASS" if report_result.ok else "FAIL",
            "abstract_words": report_result.abstract_words,
            "errors": list(report_result.errors),
            "warnings": list(report_result.warnings),
        },
    }
    console.print_json(json.dumps(payload))
    if not csv_result.ok or not report_result.ok:
        raise typer.Exit(code=2)


def _git_output(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


@app.command("submission-package")
def submission_package_command(
    csv_path: Annotated[
        Path, typer.Option("--csv", exists=True, dir_okay=False, readable=True)
    ] = Path("results/private/submission/final/jakeelamb_jl_wgs_v1.csv"),
    report_path: Annotated[
        Path, typer.Option("--report", exists=True, dir_okay=False, readable=True)
    ] = Path("results/private/submission/final/jakeelamb_track1_report.md"),
    pixi_lock: Annotated[
        Path, typer.Option("--pixi-lock", exists=True, dir_okay=False, readable=True)
    ] = Path("pixi.lock"),
    output_root: Annotated[Path, typer.Option("--output-root", file_okay=False)] = Path(
        "results/private/submission/packages"
    ),
    repository_url: Annotated[str, typer.Option("--repository-url")] = (
        "https://github.com/Jakeelamb/rare_disease"
    ),
    space_revision: Annotated[str, typer.Option("--space-revision")] = (
        "1c761cc23d90aebe6a011fd5b0b99517df42408c"
    ),
    dataset_revision: Annotated[str, typer.Option("--dataset-revision")] = (
        "f534cb0c1a607110c6dad0194299bd3dd62df542"
    ),
    allow_dirty: Annotated[bool, typer.Option("--allow-dirty")] = False,
) -> None:
    """Create an ignored content-addressed bundle from reviewed artifacts."""
    dirty = _git_output("status", "--porcelain")
    if dirty and not allow_dirty:
        console.print("[red]ERROR[/red] Git worktree is dirty; commit public-safe changes first")
        raise typer.Exit(code=2)
    package = create_submission_package(
        csv_path=csv_path,
        report_path=report_path,
        pixi_lock_path=pixi_lock,
        output_root=output_root,
        repository_url=repository_url,
        repository_commit=_git_output("rev-parse", "HEAD"),
        space_revision=space_revision,
        dataset_revision=dataset_revision,
    )
    console.print_json(
        json.dumps(
            {
                "status": "PASS",
                "bundle": str(package.bundle_dir),
                "manifest": str(package.manifest_path),
                "checksums": str(package.checksums_path),
            }
        )
    )


@app.command("render-review")
def render_review_command(
    ranked_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    output: Annotated[Path, typer.Option("--output", "-o")],
    limit: Annotated[int, typer.Option("--limit", min=1, max=10)] = 10,
) -> None:
    """Render a private, deterministic evidence card for genetics review."""
    ranked = RankedCase.model_validate_json(ranked_json.read_text(encoding="utf-8"))
    write_ranked_case_report(ranked, output, limit=limit)
    console.print(json.dumps({"report": str(output), "rendered_hypotheses": limit}, indent=2))


@app.command("candidate-read-audit")
def candidate_read_audit_command(
    ranked_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    bam: Annotated[Path, typer.Option("--bam", exists=True, dir_okay=False, readable=True)],
    json_output: Annotated[Path, typer.Option("--json-output")],
    review_output: Annotated[Path, typer.Option("--review-output")],
    limit: Annotated[int, typer.Option("--limit", min=1, max=50)] = 10,
    min_mapping_quality: Annotated[int, typer.Option("--min-mapping-quality", min=0)] = 20,
    min_base_quality: Annotated[int, typer.Option("--min-base-quality", min=0)] = 20,
    min_phase_support_fragments: Annotated[
        int, typer.Option("--min-phase-support-fragments", min=1)
    ] = 2,
) -> None:
    """Inspect leading variants directly in a local indexed BAM."""
    ranked = RankedCase.model_validate_json(ranked_json.read_text(encoding="utf-8"))
    audit = inspect_ranked_bam(
        ranked,
        bam,
        limit=limit,
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_phase_support_fragments=min_phase_support_fragments,
    )
    write_candidate_read_audit(audit, json_output=json_output, review_output=review_output)
    console.print(
        json.dumps(
            {
                "audited_alleles": len(audit.alleles),
                "audited_pairs": len(audit.pairs),
                "review": str(review_output),
            },
            indent=2,
        )
    )


@app.command("query-candidate-clinvar")
def query_candidate_clinvar_command(
    ranked_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    archive: Annotated[Path, typer.Option("--archive", exists=True, dir_okay=False, readable=True)],
    release: Annotated[str, typer.Option("--release")],
    output: Annotated[Path, typer.Option("--output", "-o")],
    hypothesis_limit: Annotated[int, typer.Option("--hypothesis-limit", min=1, max=10)] = 3,
) -> None:
    """Query a local ClinVar archive for exact leading-candidate GRCh38 alleles."""
    ranked = RankedCase.model_validate_json(ranked_json.read_text(encoding="utf-8"))
    variants = leading_variant_keys(ranked, hypothesis_limit=hypothesis_limit)
    query = query_clinvar_archive(variants, archive, release=release)
    write_clinvar_query(query, output)
    console.print(
        json.dumps(
            {
                "queried_variants": len(query.queried_variants),
                "exact_matches": len(query.matches),
                "output": str(output),
            },
            indent=2,
        )
    )


@app.command("adjudicate-leading")
def adjudicate_leading_command(
    focused_json: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    rankings: Annotated[list[str], typer.Option("--ranking", help="Repeat LABEL=RANKED_JSON")],
    recall_vcf: Annotated[
        Path, typer.Option("--recall-vcf", exists=True, dir_okay=False, readable=True)
    ],
    read_audit_json: Annotated[
        Path, typer.Option("--read-audit", exists=True, dir_okay=False, readable=True)
    ],
    clinvar_json: Annotated[
        Path, typer.Option("--clinvar", exists=True, dir_okay=False, readable=True)
    ],
    json_output: Annotated[Path, typer.Option("--json-output")],
    review_output: Annotated[Path, typer.Option("--review-output")],
    recall_label: Annotated[str, typer.Option("--recall-label")] = "DeepVariant+WhatsHap",
) -> None:
    """Integrate exact local evidence into a private manual review card."""
    focused = RankedCase.model_validate_json(focused_json.read_text(encoding="utf-8"))
    lanes: dict[str, RankedCase] = {}
    for item in rankings:
        label, separator, raw_path = item.partition("=")
        path = Path(raw_path)
        if not separator or not label or not path.is_file():
            raise typer.BadParameter(f"ranking must be LABEL=existing.json: {item}")
        if label in lanes:
            raise typer.BadParameter(f"duplicate ranking label: {label}")
        lanes[label] = RankedCase.model_validate_json(path.read_text(encoding="utf-8"))
    if not focused.candidates:
        raise typer.BadParameter("focused ranking has no candidate")
    keys = tuple(variant.key for variant in focused.candidates[0].variants)
    recall_calls = inspect_exact_vcf_calls(recall_vcf, keys, caller=recall_label)
    read_audit = CandidateReadAudit.model_validate_json(read_audit_json.read_text(encoding="utf-8"))
    clinvar = ClinVarCandidateQuery.model_validate_json(clinvar_json.read_text(encoding="utf-8"))
    adjudication = build_leading_adjudication(
        focused,
        lanes,
        recall_calls,
        read_audit,
        clinvar,
    )
    write_leading_adjudication(
        adjudication,
        json_output=json_output,
        review_output=review_output,
    )
    console.print(
        json.dumps(
            {
                "gene": adjudication.leading_candidate.gene,
                "recall_status": adjudication.recall_status,
                "phase": adjudication.recall_pair_phase,
                "review": str(review_output),
            },
            indent=2,
        )
    )
