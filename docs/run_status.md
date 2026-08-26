# Local run status

**Evidence time:** 2026-08-26 (US/Pacific)<br>
**Scope:** aggregate, non-row-level facts safe for repository review

## Current Pixi workflow

| Check | Status | Evidence |
|---|---|---|
| Locked environments | pass | `pixi install --all --locked` installed `default`, `vep`, `wgs`, `benchmark`, and `phasing` |
| Python quality | pass | 59 tests; Ruff and formatting clean; strict mypy clean across 19 source files |
| Synthetic workflow | pass | the real Snakemake `intake` and `phenotype_extraction` targets ran on generated DOCX/VCF/paired FASTQ fixtures |
| Benchmark executable path | pass | RTG Tools 3.13 vcfeval synthetic result: TP 1, FP 1, FN 0 |
| WGS command path | pass | BWA-MEM2 plus Samtools 1.23.1 aligned, name-merged, fixmated, coordinate-sorted, duplicate-marked, indexed, and quick-checked 8 synthetic reads |
| Phasing executable path | pass | WhatsHap 2.8 starts from the isolated phasing environment |
| Phase scatter/gather | pass | a synthetic five-record callset retained all records, phased four primary-chromosome variants, and preserved one non-primary variant unphased |
| Dataset payload | pass | 11 restricted artifacts; 84,985,948,316 bytes; excluded repository metadata absent |
| Real intake | pass | all compressed streams, 4 fastp lane reports, aggregate VCF inspection, and MultiQC completed in 33 minutes 18 seconds |
| Paired-read QC | review | 1,071,835,104 reads before filtering; 1,048,701,998 after filtering (97.84% retained) |
| Base quality | pass | weighted Q30 increased from 90.35% to 91.18%; 154,520,235,370 bases remained after filtering |
| Library profile | review | fastp duplication estimate 7.49-9.15%; insert-size peak 261 bp in every lane |
| Provided VCF | pass | one sample; 5,012,204 records; 4,740,790 PASS; 2,580 contigs |
| Patient-scale alignment | pass | 1,076,022,939 alignment records; 99.61% mapped; 98.12% properly paired; indexed 64.1 GB duplicate-marked BAM passed `samtools quickcheck` |
| Independent small-variant recall | pass | digest-pinned, network-disabled DeepVariant 1.10 CPU run processed 2,200,938 neural examples; 7,487,299 caller records became 7,577,588 normalized records after multiallelic decomposition |
| Whole-callset phase/gather | pass with limitation | all 7,577,588 normalized records and every per-contig count were preserved; 2,454,396 of 3,383,412 heterozygous calls (72.5%) entered 430,861 read-backed blocks |
| Leading-pair technical recall | pass | both exact private alleles are PASS heterozygous calls in supplied and DeepVariant callsets and have balanced direct BAM support |
| Leading-pair phase | unresolved | neither recall genotype has a phase set, and the candidate BAM audit found no spanning fragments; trans configuration is not established |
| Caller comparison | pass | 4,750,666 normalized alleles were shared, 56,373 were supplied-only, and 2,826,922 were DeepVariant-only; the ranking union retained disagreements for review |
| Downstream adjudication | pass | refreshed supplied and DeepVariant annotation/ranking, top-30 BAM audit, exact local ClinVar query, and final private evidence card completed |
| Submission automation | absent | no workflow rule or command creates or uploads a submission |

These checks establish execution, record preservation, and technical support;
they do not prove causality or diagnostic accuracy. The two callers use the same
library and reference, so agreement is technical replication rather than
independent biological confirmation. Short-read phase remains unresolved for
the leading private pair, and family segregation or long-read/orthogonal testing
would be needed to establish trans configuration. The current workflow ranks
SNVs and small indels and does not exclude structural, repeat, mitochondrial,
regulatory, or low-level mosaic explanations.
