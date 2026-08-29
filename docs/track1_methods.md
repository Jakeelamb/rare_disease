# Track 1 methods

**Version:** 2026-08-29
**Scope:** case-neutral public methods; patient-level inputs, evidence, and reports stay ignored

## Objective and decision unit

The workflow ranks an exact GRCh38 small variant or same-gene variant pair for
manual Track 1 review. A compound-heterozygous hypothesis is one indivisible
row: the two normalized alleles are never submitted or evaluated as unrelated
single variants. The implementation mirrors the pinned challenge scorer,
retains unique descending confidence values, and does not create or upload a
submission.

## Reproducible execution

Pixi locks five purpose-specific Linux environments (`default`, `vep`, `wgs`,
`benchmark`, and `phasing`). Snakemake 9.25.2 provides file-level dependencies,
parallel execution, and reruns after an input or command changes. Each stage is
independently runnable and ends in a private artifact that can be inspected
before continuing:

```bash
pixi run intake
pixi run phenotype-extraction
pixi run phenotype-scoring
pixi run small-variant
pixi run wgs-recall
pixi run wgs-comparison
```

The final analysis target is a review card and JSON evidence record, not a CSV.
All patient-capable paths are ignored, and `pixi run privacy-audit` checks the
Git candidate tree for genomes, alignments, phenotype documents, private
reports, credentials, and common derived formats.

## Phenotype lane

DOCX extraction is local and lossless at the source-row level. Every row starts
unresolved, and a researcher explicitly marks it proband present/absent or
family present/absent. Only reviewed proband-positive HPO terms enter scoring;
family history remains visible but cannot silently become a proband feature.

Phenotype-to-gene similarity uses the pinned HPO `v2026-06-23` ontology and
gene annotations. A disease-agnostic HPO lane is kept separate from an explicit
MVA mechanism-prior lane covering established causal genes. The mechanism prior
adds evidence but never filters the candidate universe. Agreement between the
two lanes is useful; disagreement remains a review item rather than a vote.

## Two small-variant caller architectures

The supplied VCF header identifies Sentieon Haplotyper 202308.02/GVCFtyper,
followed by GATK 4.2.4.0 hard filtering. This local-assembly/Pair-HMM-style
callset is normalized and decomposed with BCFtools before ranking.

The recall lane independently aligns all FASTQ lanes with BWA-MEM2 2.2.1 to a
source-equivalent GRCh38 no-ALT-plus-hs38d1 reference, name-merges lanes,
fixes mate tags, coordinate-sorts, marks duplicates, indexes, and quick-checks
the BAM with Samtools 1.23.1. DeepVariant 1.10.0 then applies its WGS neural
model in a digest-pinned, network-disabled CPU container. The reference and BAM
mount read-only; only the result directory is writable. Recall calls are
normalized and indexed with BCFtools 1.23.1.

Sentieon and DeepVariant use different inference machinery, so discrepancies
are informative. They still share the same biological library and reference;
caller concordance is technical replication, not independent biological
validation. The workflow starts from the union and never uses majority voting
as an inclusion rule.

## Local annotation and inheritance-aware ranking

Protein-coding exons plus 20-bp splice flanks define the fast agnostic
small-variant screen. Ensembl VEP 116.1 runs offline against its pinned GRCh38
cache and reference FASTA, retains every transcript consequence, and identifies MANE
and canonical transcripts explicitly. Native SIFT and PolyPhen values are
retained. The complete AlphaMissense GRCh38 catalogue is downloaded generically
and queried locally through the official VEP plugin; no candidate coordinate is
sent to a hosted service.

The ranker keeps heterozygous same-gene pairs and homozygous alternatives,
enforces quality and frequency thresholds when those fields are observed, and
gives missing evidence zero positive credit plus an explicit caution. Each
allele must meet a minimum consequence, clinical, or computational-support
policy. Pair-level consequence, clinical, computational, and support terms use
the weaker allele. AlphaMissense is one neural evidence type, while SIFT and
PolyPhen are averaged into one consensus so correlated predictors are not
counted as independent votes. Every score is decomposed into raw value, weight,
points, and rationale. The reported EPCR is an ordinal sorting heuristic, not a
calibrated probability.

## Phase and direct read review

Physical phase is accepted only when the VCF or read analysis records its
method. Same-phase-set `0|1`/`1|0` alleles support trans; `0|1`/`0|1` supports
cis and excludes a recessive pair. A phased genotype with no evidence method
remains unresolved. WhatsHap 2.8 phases normalized recall calls against the
indexed BAM independently by primary chromosome, then BCFtools concatenates
the shards in reference order and preserves non-primary calls unphased. The
workflow records whole-callset block statistics. Per-shard memory declarations
bound concurrency without changing chromosome-local phase semantics.

The leading focused hypotheses and top 30 agnostic hypotheses are also queried
directly in the BAM. The audit excludes unmapped, secondary, supplementary,
duplicate, QC-failed, low-mapping-quality, and low-base-quality observations.
For each exact allele it reports ref/alt/other depth, allele balance, strand,
proper-pair and soft-clipping counts, mean mapping and base quality, and mean
distance from a read end. For each pair it reports alt/alt, alt/ref, ref/alt,
ref/ref, and other shared-fragment counts. No read name or sequence is written.
An allele representation with no exact alternate read cannot be used to infer
pair phase. A direct cis/trans label requires at least two concordant
informative fragments; a lone spanning fragment remains visible as insufficient
phase support.

## Exact public-evidence adjudication

The full archived August 2026 ClinVar variant summary is downloaded generically
and streamed locally. Matches require exact assembly, chromosome, position,
reference, and alternate allele identity. The review card retains clinical
significance, review status, submitter count, evaluation date, Variation ID,
RCV accessions, and condition names; explicit no-match results are preserved.
A different nucleotide is not allowed to transfer an assertion merely because
it encodes the same amino-acid substitution.

The final private card joins:

- exact focused, agnostic, and recall ranks;
- MANE transcript, HGVSc, HGVSp, and consequence;
- supplied and recall genotypes, depth, allele depths, quality, and filters;
- direct BAM evidence and pair phase;
- population frequency, exact ClinVar assertions, and missense predictors;
- phenotype-gene and explicit mechanism evidence; and
- unresolved cautions and manual interpretation boundaries.

No automatic ACMG/AMP classification is assigned. In particular, formal PVS1
strength still requires transcript, NMD, disease-mechanism, and rescue review;
missing population frequency is not proof of absence; and computational
predictors do not establish pathogenicity.

## AI assistance and data handling

OpenAI Codex, authenticated through a ChatGPT Pro plan, assisted with coding,
testing, log review, public research, and report drafting. Bulk FASTQ, BAM, and
genome-wide VCF processing ran through local commands rather than manual web
uploads. The account's **Improve the model for everyone** setting was OFF
during the analysis, as confirmed by the entrant, so inputs and outputs were
opted out of model improvement. The case-specific submission report records
the exact assistance and data-handling boundary required by the current Track
1 instructions.

## Validation and limitations

Synthetic tests exercise phenotype curation, exact variant identity, phased
cis/trans exclusion, weakest-allele pair scoring, VEP transcript retention,
AlphaMissense parsing, direct BAM evidence, exact ClinVar matching, caller
comparison, official scorer parity, privacy, and the real Snakemake seams.
Executable smoke tests cover BWA-MEM2/Samtools alignment processing and RTG
vcfeval. Public truth-set benchmarking remains separate from case analysis.

Short-read SNV/indel analysis cannot by itself exclude structural variants,
repeat expansions, mobile elements, mtDNA variation, mosaicism below caller
sensitivity, regulatory mechanisms, or difficult mapping regions. Singleton
read phasing may leave distant pairs unresolved; family segregation or an
orthogonal assay would then be required for biological confirmation. These are
research competition predictions, not clinical diagnoses.

## Primary references

- Poplin et al., DeepVariant: <https://doi.org/10.1038/nbt.4235>
- Li, BWA-MEM2: <https://github.com/bwa-mem2/bwa-mem2>
- Martin et al., WhatsHap: <https://doi.org/10.1186/s13059-016-0918-9>
- McLaren et al., Ensembl VEP: <https://doi.org/10.1186/s13059-016-0974-4>
- Köhler et al., Human Phenotype Ontology: <https://doi.org/10.1093/nar/gkad1005>
- Cheng et al., AlphaMissense: <https://doi.org/10.1126/science.adg7492>
- Landrum et al., ClinVar: <https://doi.org/10.1093/nar/gkz972>
- Hanks et al., BUB1B and MVA: <https://doi.org/10.1038/ng1449>
- Suijkerbuijk et al., BUB1B MVA mechanism: <https://pmc.ncbi.nlm.nih.gov/articles/PMC2887387/>
