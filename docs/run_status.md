# Local run status

**Evidence time:** 2026-08-25 (US/Pacific)<br>
**Scope:** aggregate, non-row-level facts safe for repository review

## Current Pixi workflow

| Check | Status | Evidence |
|---|---|---|
| Locked environments | pass | `pixi install --all --locked` installed `default`, `vep`, `wgs`, `benchmark`, and `phasing` |
| Python quality | pass | 41 tests; Ruff clean after formatting; strict mypy clean across 17 source files |
| Synthetic workflow | pass | the real Snakemake `intake` and `phenotype_extraction` targets ran on generated DOCX/VCF/paired FASTQ fixtures |
| Benchmark executable path | pass | RTG Tools 3.13 vcfeval synthetic result: TP 1, FP 1, FN 0 |
| WGS command path | pass | BWA-MEM2 plus Samtools 1.23.1 aligned, name-merged, fixmated, coordinate-sorted, duplicate-marked, indexed, and quick-checked 8 synthetic reads |
| Phasing executable path | pass | WhatsHap 2.8 starts from the isolated phasing environment |
| Real-input DAG planning | pass | intake resolved to 8 jobs and completed; WGS recall resolves to 12 jobs in dry-run only |
| Dataset payload | pass | 11 restricted artifacts; 84,985,948,316 bytes; excluded repository metadata absent |
| Real intake | pass | all compressed streams, 4 fastp lane reports, aggregate VCF inspection, and MultiQC completed in 33 minutes 18 seconds |
| Paired-read QC | review | 1,071,835,104 reads before filtering; 1,048,701,998 after filtering (97.84% retained) |
| Base quality | pass | weighted Q30 increased from 90.35% to 91.18%; 154,520,235,370 bases remained after filtering |
| Library profile | review | fastp duplication estimate 7.49-9.15%; insert-size peak 261 bp in every lane |
| Provided VCF | pass | one sample; 5,012,204 records; 4,740,790 PASS; 2,580 contigs |
| Patient-scale WGS recall | not run in current workflow | no patient alignment or DeepVariant job was launched during the Pixi migration |
| Submission automation | absent | no workflow rule or command creates or uploads a submission |

These checks prove that the software and lightweight DAG execute. They do not
prove diagnostic accuracy.

## Historical local analysis snapshot

The following aggregate observations came from the earlier local prototype.
Its private artifacts remain ignored, but these results have **not** yet been
regenerated after the scoring, phenotype-curation, transcript, and workflow
changes in this repository.

| Analysis | Historical observation |
|---|---|
| Reference comparison | 2,580 contigs; zero header name/length differences; 5,012,204 REF alleles checked with zero mismatches |
| Prior coding/splice scan | 4,740,790 PASS calls reduced to 177,807 interval candidates; 29,351 gene-level evidence records across 10,753 genes produced 282 same-gene pair hypotheses under the old policy |
| Phase evidence | trans configuration remained unresolved in the earlier rapid and agnostic rankings |
| Orthogonal recall | reference index existed, but the full patient WGS caller lane was not run |

Treat every historical rank/count after annotation as stale until the lean Pixi
workflow reproduces it. The current ranking deliberately gives missing evidence
zero positive credit, scores pairs through the weaker allele, preserves all
transcript consequences, and does not infer a phasing method from a phased GT
alone.
