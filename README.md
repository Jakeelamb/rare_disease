# MVA Hackathon 2026: transparent rare-disease research

This is a Python/Snakemake workspace for the Sage Bionetworks **Rare Disease,
Real Kid: MVA Hackathon 2026**. It is built for fast, inspectable research:
prototype a lane, compare it with the others, correct mistakes, and rerun.

The repository does **not** generate or upload a submission. Patient data and
all row-level results stay under ignored local directories.

## What is implemented

- one locked Pixi workspace with isolated default, VEP, WGS, benchmark, and
  phasing environments;
- local FASTQ/VCF QC and phenotype extraction;
- explicit manual phenotype curation before HPO-to-gene scoring;
- a disease-prior lane and an agnostic coding/splice lane;
- typed variant evidence, exact same-gene pair construction, visible score
  decomposition, missing-evidence handling, and phase cautions;
- union-first comparisons across ranking lanes and advisory truth-set metrics
  across callers;
- a repaired BWA-MEM2/Samtools/DeepVariant/WhatsHap WGS lane;
- synthetic workflow, RTG vcfeval, WGS-command, and phasing smoke tests.

The Pixi migration is tested on synthetic data, and the real intake/QC stage
has completed against the full local challenge payload. The full patient-scale
WGS recall has not been run through this version of the workflow. Current and
historical aggregate evidence are separated in [run status](docs/run_status.md).

## Pixi setup

Install Pixi, make it available on `PATH`, then materialize every locked
environment:

```bash
pixi --version
pixi install --all --locked
```

`pixi.lock` is the Linux x86-64 dependency contract. Do not maintain a second
environment definition beside it.

Run the complete fast acceptance suite with:

```bash
pixi run check
pixi run -e benchmark benchmark-tool-smoke
pixi run -e wgs wgs-tool-smoke
pixi run -e phasing phasing-tool-smoke
```

## Restricted dataset

Authentication stays outside the repository. Never paste a token into a
command, notebook, config file, issue, or commit.

```bash
pixi run hf-login
pixi run hf-whoami
pixi run download-data
```

`download-data` is pinned to dataset revision
`f534cb0c1a607110c6dad0194299bd3dd62df542`, writes to ignored `data/`, and
excludes the dataset README and `.gitattributes`. Hugging Face resumes files
already present in that local directory.

Before any Git operation that could publish files, run:

```bash
pixi run privacy-audit
```

The audit fails closed on controlled-data paths, patient-capable genomic file
types, high-confidence credential patterns, and any Git candidate larger than
5 MiB. Large reference databases and all challenge inputs stay local and
ignored.

## Research stages

Outputs go to ignored `results/private/work/`. These are ordinary, rerunnable
Snakemake targets—not approvals or promotion gates. The default target is the
lightweight intake/QC stage.

```bash
pixi run intake
pixi run phenotype-extraction
```

Inspect:

- `results/private/work/qc/multiqc/multiqc_report.html`
- `results/private/work/qc/vcf_summary.json`
- `results/private/work/phenotype/phenotype.review.md`

Automatic extraction deliberately leaves every phenotype row unresolved.
Create `phenotype.curated.json` by assigning exactly one decision to every row:
`proband_present`, `proband_absent`, `family_present`, or `family_absent`.

```bash
pixi run mva curate-phenotype \
  results/private/work/phenotype/phenotype.json \
  --decision 1=proband_present \
  --decision 2=family_present \
  --reviewer "RESEARCHER NAME" \
  --output results/private/work/phenotype/phenotype.curated.json
```

Repeat `--decision` for every extracted row, then continue:

```bash
pixi run phenotype-scoring
pixi run small-variant
```

The main review artifacts are:

- `ranking/mva_prior.review.md`
- `ranking/agnostic_coding.review.md`
- `ranking/small_variant.lane_comparison.md`

all beneath `results/private/work/`. A disagreement is preserved for research
review; the workflow never resolves it by majority vote.

## Comparing algorithms

RTG Tools 3.13 is isolated in the `benchmark` environment. The generic adapter
in `scripts/run_vcfeval.py` writes standardized TP/FP/FN, runtime, and memory
metrics. Place baseline, candidate, and union metrics at the configured private
paths, then run:

```bash
pixi run benchmark
```

The resulting `benchmark/tool_comparison.json` is advisory. It reports rescued
truth variants, added false positives, and runtime ratio; it does not block any
other experiment.

The orthogonal WGS lane is also directly runnable:

```bash
pixi run wgs-smoke
pixi run wgs-recall
pixi run wgs-comparison
```

`wgs-recall` is the expensive patient-scale path. It aligns all lanes, marks
duplicates, runs the digest-pinned DeepVariant CPU container without network
access, then normalizes and read-phases calls. Run it only when that compute is
intended. `wgs-comparison` is a separate downstream target that requires the
manually curated phenotype, ranks DeepVariant coding/splice candidates, and
writes the supplied-callset-versus-DeepVariant discrepancy report.

## Manual submission boundary

No command or workflow rule creates a submission CSV or uploads anything. A
read-only schema checker remains available for a CSV that researchers prepare
manually:

```bash
pixi run mva validate-submission path/to/manual.csv
```

## Privacy and scientific scope

The input genome and phenotype belong to a real child. Do not redistribute raw
or derived patient data, send patient-level content to hosted APIs, or commit
private reports. Follow the official controlled-access, no-recontact, embargo,
and deletion requirements summarized in
[challenge contract](docs/challenge_contract.md).

This is research software, not a clinical diagnostic system. See the
[architecture](docs/architecture.md) and pinned
[challenge and methods evidence brief](docs/research/challenge_and_methods.md)
for assumptions, sources, and unresolved limitations.
