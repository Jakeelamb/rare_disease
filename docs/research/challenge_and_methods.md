# Rare Disease, Real Kid: MVA Hackathon 2026

## Challenge evidence and method brief

**Evidence snapshot:** 2026-08-25<br>
**Official Space revision:** `37e25dceda63ecec7c5b2ebeffd1ea0052ad886e`<br>
**Official dataset revision:** `f534cb0c1a607110c6dad0194299bd3dd62df542`

This note is pinned to the revisions above because the challenge launched on 2026-08-24 and its rules, code, and timeline may still change. Re-check the live Space before every submission.

## Executive conclusions

1. This is a single-proband, two-track challenge built around a real child with Mosaic Variegated Aneuploidy (MVA). Track 1 asks for a ranked causal-variant list; Track 2 asks for mechanism-grounded hypotheses about already approved medicines. Track 2 proposals are explicitly hypotheses for follow-up, not evidence of efficacy or medical advice. [Official overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/about.py)
2. The public Track 1 evaluator says the clinical answer is a **compound-heterozygous pair**. Full rank credit requires putting both exact GRCh38 variants in the same CSV row. [Evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)
3. The provided VCF is a useful starting artifact, not a complete diagnostic workflow: its header identifies Sentieon Haplotyper/GVCFtyper followed by GATK VariantFiltration, and it contains a single sample with standard small-variant genotype/quality fields. There is no functional annotation in the observed header and no provided SV/CNV callset.
4. The official repositories contain the submission UI, scorer, and templates, but no starter diagnostic notebook, baseline analysis pipeline, trained model, or benchmark output. A transparent workflow therefore needs to establish its own provenance, annotation, prioritization, review, and exact-output validation.
5. The safest high-performing design is a **two-lane analysis**: an agnostic phenotype/inheritance-aware whole-genome ranking plus an explicit MVA-mechanism hypothesis lane. The latter must not hard-filter away novel or noncanonical explanations.
6. Do not send patient-level clinical or variant data to hosted LLM/API providers until the organizers answer the open privacy question. The published rules prohibit granting data access and require deletion of source and derived data; organizers had not answered the specific third-party API question as of this snapshot. [Open discussion #2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)

## 1. Official objective, tracks, and schedule

### Track 1: Variant Prediction

Goal: use the proband's FASTQs, VCF, and clinical phenotype to identify the specific variant or variants driving the condition. Submit a ranked list, public GitHub repository, and methods report. Secondary/incidental findings may be reported but are qualitatively reviewed. The automated result is displayed on a public leaderboard. [Overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/about.py) · [Track 1 instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track1.py)

### Track 2: Drug Repurposing

Goal: use the genetic result and its mechanism to identify already market-approved medicines whose mechanism may plausibly target the disrupted pathway. The report must characterize loss/gain of function, the affected pathway, and downstream consequences. Submit one report (PDF or Markdown), a public GitHub repository, and a three-minute YouTube/Vimeo pitch. [Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py) · [Track 2 instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track2.py)

### Published schedule

| Milestone | Published date |
|---|---:|
| Dataset available | 2026-08-24 |
| Space says submissions open | 2026-08-25 |
| Submission close / Track 1 freeze | 2026-10-24 23:59 UTC |
| Published judging window | 2026-10-24 through 2026-11-24 |
| Published winner announcement | 2026-11-25 |

The Space says dates may change. The dataset README instead describes the submission period as 2026-08-24 through 2026-10-24. [Overview timeline](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/about.py) · [Dataset README](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/blob/f534cb0c1a607110c6dad0194299bd3dd62df542/README.md)

### Prize language

The Space lists a $50,000 total pool: $25,000 cash plus $25,000 Claude credits, divided into first ($12k + $12k), second ($7k + $7k), third ($4k + $4k), and innovation/community ($2k + $2k). It does not clearly publish whether these are overall awards or how the two tracks map to prizes. [Official rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py)

## 2. Dataset inventory and schemas

The public tree at the pinned dataset revision contains 11 challenge artifacts plus the README and `.gitattributes`; total repository payload is 84,985,955,953 bytes (approximately 85.0 GB decimal). [Pinned dataset tree](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/tree/f534cb0c1a607110c6dad0194299bd3dd62df542)

| Artifact | Count | Exact bytes | Observed role/schema |
|---|---:|---:|---|
| `Challenge_Clinical_Phenotype_1.docx` | 1 | 16,865 | Confidential participant reference. One clinical-feature table with columns `Clinical Feature`, `HPO Term`, `HPO ID`, and `Presentation / Notes`, plus interpretive notes. Do not duplicate its row-level content into a public repository. |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | 1 | 315,153,971 | bgzip VCF 4.2, GRCh38, one sample (`WGS_EX2312012`), contigs named `1`...`22`,`X`,`Y`,`M` (no `chr` prefix), standard `GT`, `AD`, `DP`, `GQ`, `PGT`, `PID`, `PL` fields. Header records Sentieon Haplotyper 202308.02, Sentieon GVCFtyper, and GATK 4.2.4.0 hard filtering. |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` | 1 | 2,343,376 | Tabix index for the VCF. |
| `..._L001_...fastq.gz` through `..._L004_...fastq.gz` | 8 | 84,668,434,104 total | Four lanes of paired-end reads (`R1`/`R2` for lanes `L001`-`L004`). Individual files are 10.20-11.04 GB. |

The gated phenotype document inspected locally has SHA-256 `0b8129496f239d766c09506f86c23ca2392d7f5f1db13ee8ee513688fea584b7`; the VCF index has SHA-256 `6f8fed62f11c475fc63a8e2b50925ffc7be33b6930225e821cb977951761a2e0`. These hashes identify the inputs without redistributing clinical contents.

### Important data facts

- The dataset is one subject, not a train/test cohort. There are no provided training labels.
- The public scorer describes the hidden clinically validated ground truth as two causal variants in a compound-heterozygous configuration. It does not disclose their coordinates. [Evaluator header](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)
- The distributed VCF is a small-variant callset. The raw FASTQs permit an orthogonal re-call and SV/CNV analysis, but neither a BAM/CRAM nor an SV/CNV VCF is present at the pinned revision.
- The rules describe phenotype data as standardized HPO terms and data as VCF with optional BAM/CRAM, but the actual release is a DOCX HPO table, VCF/TBI, and FASTQs. The repository tree is authoritative for implementation. [Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py)

## 3. Track 1 submission contract and exact scorer behavior

### CSV schema

The required template has 12 columns:

```text
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

Constraints:

- `proband_id` must be `PROBAND01`.
- Coordinates must be GRCh38.
- One row represents either one candidate variant or one proposed compound-heterozygous pair.
- `epcr` is an estimated probability of causal relationship in `(0,1]`; the scorer sorts rows descending by this value.
- `finding_type` must be `primary` or `secondary`.
- At most 10 rows are accepted.
- Each submission also requires a public GitHub URL and a PDF/Markdown report.
- Each authenticated Hugging Face participant receives six Track 1 uploads; only the highest-scoring team entry is shown.

[Template](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/static/templates/track1_submission_template.csv) · [Upload validation](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track1.py)

### Automated metrics

For a full exact pair match, rank points are:

| Sorted rank | Points |
|---:|---:|
| 1 | 100 |
| 2-3 | 50 |
| 4-5 | 25 |
| 6-10 | 10 |

If no full pair is present but a row overlaps one of the two true variants, rank points are half the corresponding tier. F-max is the maximum variant-level F1 score across all unique EPCR thresholds, where all variants in rows at or above a threshold form the predicted set. [Exact scorer](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)

The metric is adapted from the CAGI6 Rare Genomes Project assessment, but the hackathon code is the binding implementation. The source paper is Stenton et al., 2024, [“Critical assessment of variant prioritization methods for rare disease diagnosis within the rare genomes project”](https://doi.org/10.1186/s40246-024-00604-w).

### Acceptance-test implications

These are direct consequences of the code, not speculative tactics:

1. **Submit the two candidate alleles together.** Two separate single-variant rows can recover both alleles for F-max, but cannot earn full-match rank points.
2. **Normalize exactly.** The scorer uppercases REF/ALT but does not left-normalize, trim alleles, reconcile multiallelic records, or normalize chromosome names. The provided VCF uses `1`, whereas the template and scorer fallback use `chr1`; submission generation must intentionally emit the challenge's `chrN` convention and validate it locally against a synthetic key.
3. **Put the strongest pair first.** A correct pair at rank 1 earns 100 rank points and can reach F-max 1.0 if false variants are below its threshold. False variants above the true pair lower both rank and achievable F-max.
4. **`finding_type` is not a scoring exclusion.** Despite prose saying secondary findings will not hurt automated scoring, the evaluator ignores `finding_type` when sorting and computing F-max. A high-EPCR secondary row can push the causal pair down and enter the predicted set. Keep incidental findings below primary candidates or separate them from the scored CSV pending organizer clarification.
5. **No combined leaderboard formula is published.** The leaderboard exposes Rank Points and F-max separately; the final tie-break, weighting, and role of qualitative method review are not specified.

## 4. Track 2 submission and judging contract

Required artifacts:

- one PDF or Markdown report;
- one public GitHub repository;
- one three-minute YouTube or Vimeo pitch;
- optional judge notes.

The panel weights are:

| Criterion | Weight | Published question |
|---|---:|---|
| Scientific rigor | 35% | Is mechanism characterization sound and is the candidate supported by it? |
| Potential impact | 25% | Could validation advance understanding or diagnosis for this child or others with MVA? |
| Innovation | 25% | Is the angle, method, or tool genuinely creative? |
| Scalability | 15% | Could the approach generalize beyond this case? |

[Overview and weights](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/about.py) · [Track 2 upload contract](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track2.py)

The UI says one submission per team and asks one designated member to upload. The backend enforces one per Hugging Face username, while the prose says duplicate submissions from teammates will not be reviewed. Treat the team limit as binding.

## 5. Rules, privacy, licensing, and publication restrictions

### Binding participant obligations

- Participants must be at least 18 and every team member must register and accept the rules.
- No recontact of the data subject, family, or MVA Society contacts.
- No release or access grant to anyone; every team member must be registered.
- Safeguards are required, and suspected unauthorized disclosure must be reported to Sage's Privacy and Compliance Office.
- All source data and **intermediate or derived datasets** must be deleted from local systems, cloud, notebooks, and private repositories within 30 days after hackathon close.
- Participants must email `RarediseaserealkidMVAhackathon2026@synapse.org` to confirm deletion.
- Participant submissions may be rerun and are released under CC BY with participant attribution.
- A publication embargo begins at challenge close and ends when the organizers publicly post their summary report or preprint. During it, participants cannot submit a peer-reviewed manuscript based on the dataset; preliminary abstracts/posters require prior written approval.
- Public communications must use the required acknowledgement and must not add re-identifying information.

[Official rules at pinned revision](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py)

### License interpretation

The dataset card declares `cc-by-4.0`, but the same card gates access behind special attestations, and the rules/FAQ prohibit redistribution and require deletion. Participant code, reports, predictions, and results are expressly CC BY 4.0; patient-level source data remain controlled by the special hackathon terms. Until organizers resolve the metadata tension, the conservative rule is: **never commit source data or row-level extracts, and do not treat the dataset-card license as permission to redistribute patient data.** [Dataset card](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/blob/f534cb0c1a607110c6dad0194299bd3dd62df542/README.md) · [FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/faq.py)

### Local/hosted model boundary

An organizer answer is still pending on whether variant rows, HPO terms, candidate lists, or prompts may be sent to hosted LLM APIs and whether provider logs fall under the deletion obligation. Until answered, use local inference only for patient-derived content, or restrict hosted tools to public literature and non-patient-specific questions. [Open organizer discussion](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)

## 6. What the official baseline does and does not provide

At the pinned revision, the official Space is a Gradio application with only `gradio`, `huggingface_hub`, and `python-dotenv` as declared dependencies. Its public technical assets are:

- submission templates and method-description workbook;
- exact Track 1 CSV parser and evaluator;
- ground-truth loader whose production key is in a private dataset;
- submission storage and leaderboard UI;
- rules, FAQ, and Track 2 upload UI.

There is no official diagnostic baseline code, environment, notebook, model, candidate ranking, expected runtime, or expected cost. No other SageBio model/dataset/Space matching this MVA challenge was visible through the Hugging Face organization APIs on 2026-08-25. [Space repository tree](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/tree/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e) · [Requirements](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/requirements.txt)

The provided VCF is therefore the organizer's only analysis starting point. Its header provenance should be preserved and compared against, not silently replaced.

## 7. Evidence-backed Track 1 workflow

### Design principle

Use Python for orchestration, evidence tables, reports, and tests; call established version-pinned genomics tools for alignment/calling/annotation. Reimplementing mature callers in Python would reduce transparency and scientific defensibility.

### Stage A: immutable intake and QC

1. Record dataset revision, file size, SHA-256, acquisition time, and access-control status in a machine-readable manifest.
2. Validate FASTQ gzip integrity, paired read counts, read-name pairing, lane identity, base quality, adapter content, GC distribution, and duplication.
3. Parse and preserve the provided VCF header, caller versions, reference identity, sample name, filter definitions, contig convention, and record counts.
4. Keep all patient-derived data beneath a gitignored controlled root; publish only code, environment locks, schemas, aggregate runtime/QC, and appropriately minimized evidence.

### Stage B: structured phenotype model

Extract the gated DOCX table into an internal HPO manifest with positive observations, onset/context, and family-history attribution. Preserve exact source-cell provenance. Run two complementary phenotype models:

- **Exomiser** for inheritance-aware combination of variant pathogenicity, frequency, and phenotype similarity. The original validation showed that human/model-organism phenotype integration improves variant prioritization, including cases without parental sequence. [Smedley et al., 2015](https://www.nature.com/articles/gim2015137)
- **LIRICAL** for an interpretable likelihood-ratio differential with per-phenotype contributions and post-test probabilities. [Robinson et al., 2020](https://doi.org/10.1016/j.ajhg.2020.06.021)

HPO is the standard computable phenotype vocabulary used by both. [Köhler et al., 2021](https://doi.org/10.1093/nar/gkaa1043)

Keep an **agnostic lane** and an **MVA-prior lane**. MVA is known to result from biallelic disruption of mitotic/chromosome-segregation genes including `BUB1B`, `CEP57`, and `TRIP13`, but these genes are priors, not filters:

- Biallelic `BUB1B` variants were established in MVA with childhood cancer predisposition. [Hanks et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15475955/)
- Biallelic loss-of-function `CEP57` variants were established as another MVA cause. [Snape et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3508359/)
- Biallelic loss-of-function `TRIP13` variants were shown to impair the spindle assembly checkpoint and cause chromosome missegregation/cancer predisposition. [Yost et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5493194/)

### Stage C: normalize and annotate the provided callset

1. Split multiallelic sites and left-normalize against the exact GRCh38 reference used for analysis.
2. Annotate all relevant transcripts/consequences with Ensembl VEP, retaining transcript IDs and versions rather than collapsing immediately to one consequence. VEP is designed for coding and noncoding variant annotation and prioritization. [McLaren et al., 2016](https://doi.org/10.1186/s13059-016-0974-4)
3. Add versioned population frequency and constraint evidence (gnomAD), clinical assertions with review status/conflicts (ClinVar), predicted loss-of-function confidence, splice effect (SpliceAI), and missense evidence (AlphaMissense plus calibrated ensemble scores where available).
4. Never use any one computational score as a pathogenicity verdict:
   - gnomAD provides population variation and gene constraint at large scale. [Karczewski et al., 2020](https://www.nature.com/articles/s41586-020-2308-7)
   - ClinVar is an archive of submitted variant-condition interpretations and supporting evidence, not a single infallible label. [Landrum et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29165669/)
   - SpliceAI predicts splice disruption from primary sequence and experimentally validated many cryptic-splice predictions. [Jaganathan et al., 2019](https://doi.org/10.1016/j.cell.2018.12.015)
   - AlphaMissense predicts missense pathogenicity proteome-wide; its scores remain one evidence type. [Cheng et al., 2023](https://doi.org/10.1126/science.adg7492)

### Stage D: inheritance-aware pair generation

Because the public evaluator says the answer is compound heterozygous, the central ranking unit must be a **pair**, not an isolated variant.

For each gene:

1. generate all plausible heterozygous rare-variant pairs;
2. retain each allele's QC, consequence, frequency, clinical evidence, phenotype relevance, and independent read support;
3. score compatibility with autosomal-recessive disease and the expected molecular mechanism;
4. penalize pairs likely to be in cis, recurrent technical artifacts, common haplotypes, or pairs whose combined genotype is incompatible with known disease mechanism;
5. explicitly label phase as read-backed, statistical, family-backed, or unresolved.

The single-proband release has no parental genotypes, so trans configuration cannot generally be assumed. Read-backed phasing can help when variants share informative fragments; HapCUT2 is a published haplotype-assembly method across sequencing technologies, but short-read phase blocks remain limited. [Edge et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411775/)

### Stage E: orthogonal raw-read analysis

The provided Sentieon/GATK VCF should remain one independent lane. Reprocess raw reads to test caller sensitivity and uncover variant classes absent from it:

1. lane-aware QC and alignment to the exact pinned GRCh38 build;
2. an orthogonal small-variant callset such as DeepVariant, whose published model calls SNPs/small indels from aligned reads; [Poplin et al., 2018](https://www.nature.com/articles/nbt.4235)
3. SV/indel analysis from paired/read-split signals (for example Manta); [Chen et al., 2016](https://doi.org/10.1093/bioinformatics/btv710)
4. independent read-depth CNV/aneuploidy analysis suitable for a single genome (for example CNVnator), interpreted cautiously without matched controls; [Abyzov et al., 2011](https://doi.org/10.1101/gr.114876.110)
5. candidate-level pileup and alignment review for every final allele.

DeepVariant, Manta, and CNVnator address different variant classes; their outputs should not be blended into a score without class-specific QC. The confirmed scoring target is a small-variant pair, so SV/CNV work is primarily a safeguard, mechanistic context, and possible innovation evidence rather than a reason to neglect pair ranking.

### Stage F: evidence review and ranking

For every top pair, generate a transparent evidence card:

- normalized genomic and transcript notation;
- genotype quality, depth, allele balance, strand/read-position evidence, and caller concordance;
- population frequency by ancestry and quality flags;
- consequence for every relevant transcript;
- known disease/gene mechanism and inheritance;
- ClinVar assertions, review levels, conflicts, and dates;
- HPO similarity plus per-feature contribution;
- functional literature with direct quotations kept below copyright limits and linked to primary sources;
- ACMG/AMP evidence codes, strength, rationale, and reviewer;
- phase status and validation needed;
- reasons the pair could be wrong.

Use ACMG/AMP terminology and evidence categories for transparent interpretation, while clearly marking the workflow as research and the result as requiring clinical confirmation. [Richards et al., 2015](https://www.nature.com/articles/gim201530)

### Stage G: output validation

Before consuming one of six uploads:

1. validate exactly 12 required columns, `PROBAND01`, at most 10 rows, EPCR `(0,1]`, and descending unique probabilities;
2. assert `chr` chromosome convention and GRCh38 normalized REF/ALT against a pinned reference;
3. assert both alleles of each compound-heterozygous hypothesis occupy the same row;
4. run a vendored, revision-pinned copy of the official evaluator against synthetic known pairs, including indels and chromosome-prefix cases;
5. prepare the CSV manually from the reviewed evidence table, then run the
   read-only contract checker; the analysis workflow must neither serialize nor
   upload it;
6. archive the code commit, environment lock, manifest, report hash, CSV hash, and exact scorer test log.

## 8. Evidence-backed Track 2 workflow

Do not begin with a drug list. Begin with a falsifiable mechanism chain:

```text
variant pair -> transcript/protein effect -> cellular direction of effect
-> pathway/cell phenotype -> intervention direction -> approved drug
-> exposure/tissue plausibility -> pediatric safety -> validation experiment
```

### Candidate-generation lanes

1. **Direct target modulation:** approved agonist/antagonist/stabilizer/replacement acting on the causal gene product, if directionally appropriate.
2. **Pathway compensation:** modulate an upstream/downstream node expected to restore the measured defect without worsening chromosomal instability.
3. **Phenotypic signature reversal:** only if a disease-relevant expression or perturbation signature exists. Connectivity Map/L1000 supports perturbational matching, but the source release contains no transcriptome, so signature reversal would be hypothesis-generating rather than patient-specific proof. [Subramanian et al., 2017](https://doi.org/10.1016/j.cell.2017.10.049)
4. **Cancer-risk/symptom management:** keep supportive or surveillance concepts separate from disease-modifying claims.

### Evidence sources and provenance

- Use Open Targets to assemble target-disease, direction-of-effect, clinical precedence, tractability, and safety evidence with source-level provenance. [Open Targets Platform, 2025](https://doi.org/10.1093/nar/gkae1128)
- Use ChEMBL for curated target/activity/mechanism records, preserving assay type, organism, target confidence, activity units, and source document. [Zdrazil et al., 2024](https://doi.org/10.1093/nar/gkad1004)
- Confirm **market approval and current labeling** in regulator-owned sources such as [Drugs@FDA](https://www.fda.gov/drugs/drug-approvals-and-databases/about-drugsfda) and [DailyMed](https://dailymed.nlm.nih.gov/dailymed/), not from a knowledge-graph edge alone.
- Search [ClinicalTrials.gov](https://clinicaltrials.gov/) for disease/pathway experience and [PubMed](https://pubmed.ncbi.nlm.nih.gov/) for original mechanistic, cellular, animal, pharmacokinetic, and safety evidence.
- Record query date, database release/API response hash, identifiers, and every transformation so the ranking can be rerun.

### Drug evidence card

Each proposed medicine should have:

1. generic name, regulator, jurisdiction, approval status, label date, indication, formulation, and age restrictions;
2. direct molecular target(s), action type, potency range, assay context, and evidence source;
3. explicit desired direction of modulation and why it compensates the variant mechanism;
4. MVA-relevant cell/tissue plausibility and exposure constraints;
5. known adverse effects, boxed warnings, cancer/genome-instability concerns, interactions, and pediatric uncertainties;
6. evidence tier: direct human, disease-model, cellular, target-only, network-only, or computational;
7. falsification criteria and a staged validation experiment with readouts, controls, dose range, and stop conditions;
8. limitations and an explicit “not for clinical use” statement.

For a chromosome-instability disorder, any intervention affecting mitosis, DNA damage, apoptosis, or oncogenic selection requires especially strong safety analysis. A computationally attractive reversal score is not enough.

## 9. Transparency and reproducibility requirements

The public repository should contain no controlled patient data. It should contain:

- a pinned Pixi environment matrix and platform lock backed by Conda-Forge/Bioconda;
- workflow DAG and one command per documented entry point;
- immutable input-manifest schema with checksums represented by safe placeholders/example data;
- configuration separated from code;
- provenance for tool/container/database versions and reference hashes;
- unit tests for normalization, pair construction, EPCR ordering, schema, and scorer parity;
- integration tests on public synthetic/GIAB-style data;
- structured evidence tables and deterministic report generation;
- resource/runtime/cost capture for every task;
- a data-retention/deletion manifest that covers caches, temporary files, cloud objects, notebooks, logs, and derived artifacts;
- a model card/method report stating what was automated, what was manually reviewed, public versus proprietary evidence, strengths, limitations, and whether compound-heterozygous pairs are supported.

The official methods workbook explicitly asks for automation versus manual curation, public/proprietary data, compound-heterozygous capability, secondary findings, runtime/cost, and strengths/limitations. Those fields should be generated from real logs rather than estimated after the fact. [Methods workbook](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/static/templates/methods_description_form.xlsx)

## 10. Unresolved blockers and organizer questions

1. **Hosted API/LLM privacy:** no organizer answer yet on whether patient-derived HPO/variant data can be sent to third-party APIs, what retention terms are acceptable, or how provider logs interact with deletion. [Discussion #2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)
2. **Dataset license versus special terms:** the card says CC BY 4.0 while rules prohibit redistribution and mandate deletion. Ask organizers to state the governing license/data-transfer terms unambiguously.
3. **Derived-data boundary:** rules require deletion of all intermediate/derived datasets, but also allow public code, models, and derived outputs. Ask which artifacts—HPO lists, candidate tables, annotated VCFs, prompts, logs, and evidence cards—may be retained or published.
4. **Secondary-finding scorer mismatch:** prose says secondary rows do not hurt, but the evaluator includes all rows regardless of `finding_type`. Ask whether the code or documentation will change.
5. **Chromosome normalization:** provided VCF contigs omit `chr`; template and evaluator fallback use `chr`. Ask organizers to confirm the production ground-truth convention.
6. **Track 1 ranking:** no published formula combines Rank Points and F-max, no tie-break, and no precise role for qualitative method review in final awards.
7. **Prize allocation:** no explicit mapping of awards to Track 1 versus Track 2.
8. **Timeline conflict:** rules/FAQ say approximately 2-3 months of judging, while the dated timeline shows approximately one month and a 2026-11-25 announcement.
9. **Required pitch by track:** rules broadly say every submission includes a three-minute pitch, while Track 1's upload UI does not request a video and Track 2 does. Ask whether Track 1 needs a pitch elsewhere.
10. **Dataset citation:** rules require a reference from a Synapse page, but no Synapse dataset citation/link is published in the inspected Space or dataset README.
11. **Current local acquisition:** all 11 requested challenge files were present on 2026-08-25 and Hugging Face dry-run reported 0 bytes remaining. The VCF contains 5,012,204 records (4,740,790 `PASS`) and one sample. All compressed streams passed integrity validation; read-level quality profiling remains in progress.

## 11. Source index

### Official challenge sources

- [Challenge Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
- [Pinned Space repository](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/tree/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e)
- [Pinned dataset repository](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/tree/f534cb0c1a607110c6dad0194299bd3dd62df542)
- [Official rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py)
- [Exact evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)
- [Track 1 submit code/instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track1.py)
- [Track 2 submit code/instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/submit_track2.py)

### Primary method and disease sources

- Stenton et al. 2024, CAGI6 Rare Genomes Project: <https://doi.org/10.1186/s40246-024-00604-w>
- Smedley et al. 2015, Exomiser: <https://www.nature.com/articles/gim2015137>
- Robinson et al. 2020, LIRICAL: <https://doi.org/10.1016/j.ajhg.2020.06.021>
- Poplin et al. 2018, DeepVariant: <https://www.nature.com/articles/nbt.4235>
- McLaren et al. 2016, Ensembl VEP: <https://doi.org/10.1186/s13059-016-0974-4>
- Karczewski et al. 2020, gnomAD constraint: <https://www.nature.com/articles/s41586-020-2308-7>
- Landrum et al. 2018, ClinVar: <https://pubmed.ncbi.nlm.nih.gov/29165669/>
- Jaganathan et al. 2019, SpliceAI: <https://doi.org/10.1016/j.cell.2018.12.015>
- Cheng et al. 2023, AlphaMissense: <https://doi.org/10.1126/science.adg7492>
- Chen et al. 2016, Manta: <https://doi.org/10.1093/bioinformatics/btv710>
- Abyzov et al. 2011, CNVnator: <https://doi.org/10.1101/gr.114876.110>
- Edge et al. 2017, HapCUT2: <https://doi.org/10.1101/gr.213462.116>
- Richards et al. 2015, ACMG/AMP interpretation: <https://www.nature.com/articles/gim201530>
- Hanks et al. 2004, `BUB1B` and MVA: <https://pubmed.ncbi.nlm.nih.gov/15475955/>
- Snape et al. 2011, `CEP57` and MVA: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3508359/>
- Yost et al. 2017, `TRIP13` and MVA: <https://pmc.ncbi.nlm.nih.gov/articles/PMC5493194/>
- Buniello et al. 2025, Open Targets Platform: <https://doi.org/10.1093/nar/gkae1128>
- Zdrazil et al. 2024, ChEMBL: <https://doi.org/10.1093/nar/gkad1004>
- Subramanian et al. 2017, Connectivity Map/L1000: <https://doi.org/10.1016/j.cell.2017.10.049>
