# State of the art for a transparent singleton short-read WGS rare-disease workflow (2026)

**Research date:** 2026-08-25<br>
**Scope:** one proband, short-read whole-genome sequencing, no parental samples, with compound-heterozygous hypotheses central to the competition.<br>
**Evidence policy:** primary papers, official tool documentation/repositories, authoritative resource pages, and immutable official challenge source only.<br>
**Privacy boundary:** this review did **not** inspect patient phenotype rows, variant rows, candidates, genes, coordinates, alignments, or derived patient-level files. It is a workflow and methods review, not a case analysis.

## How to read this note

- **Established** means directly supported by a primary paper or official documentation.
- **Inference** means a consequence of the public competition contract or of the limitations of singleton short-read WGS.
- **Recommendation** is the proposed decision for this competition; it is not a clinical claim.

Version and license statements are a snapshot as of the research date. Tool code, model weights, annotation files, and bundled third-party resources can have different terms. “Locally usable” means technically able to run without uploading the proband after resources are downloaded; it does not imply unrestricted commercial use.

**Release-label discipline.** Throughout this note, “current release” means a tagged artifact on the official release page. A version displayed by `development`/`latest` documentation is not called released unless the official release page also contains it. In particular, Exomiser **15.1.0** is the current tagged application; Read the Docs development pages labeled 15.1.1 are forward-looking documentation, not proof of a 15.1.1 release.

## Executive decision

The most defensible winning workflow is not a maximal ensemble. It is a staged system with a narrow, auditable backbone:

1. lock the reference/contig/normalization contract and prove BAM integrity;
2. use one current high-accuracy short-variant caller as the primary callset (DeepVariant WGS is the strongest default), with GATK HaplotypeCaller as an independently implemented assembly/Bayesian rescue caller;
3. construct recessive hypotheses at the **gene and allele-pair level**, not as a ranked bag of variants;
4. combine explicit phenotype similarity, disease inheritance, population frequency, predicted molecular consequence, read evidence, and phase status into a decomposable score;
5. run read-backed phasing for candidate pairs, but represent unresolved phase honestly—short reads cannot usually prove trans for distant variants;
6. add orthogonal lanes for CNV/SV, repeat expansion, mtDNA, and mobile elements in a deliberate second stage;
7. use deep-intronic splice scoring only within phenotype-supported genes or as a partner to a strong coding allele;
8. calibrate the entire ranking and submission transformation on public truth/synthetic cases, never on the hidden patient answer.

This recommendation follows the public evaluator. Its immutable source says the answer is a clinically validated compound-heterozygous pair, accepts at most ten rows, sorts by EPCR, awards full rank credit only when the exact two-variant set occurs in one row, and computes F-max over the union of individual variants above each EPCR threshold ([official evaluator at Space commit `37e25dc`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)). The official template requires both alleles in the same row ([submission template](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/static/templates/track1_submission_template.csv)).

**Inference:** the competition unit is a *pair hypothesis*. Optimizing variant-level pathogenicity independently is structurally misaligned with the score. False partners also damage F-max, so a short, calibrated list is preferable to ten weak pair permutations.

## 1. Competition-specific objective and constraints

### Public scoring contract

**Established.** The official evaluator:

- parses variants as literal `(chrom, pos, upper(ref), upper(alt))` tuples;
- gives a full match only for set equality with the hidden two-variant truth;
- gives half-weight rank credit if there is no full pair but a row intersects the truth;
- permits at most ten rows;
- sorts by descending EPCR, using input order to break ties;
- treats `finding_type` as informational in automated scoring; and
- finds the maximum F-measure over distinct submitted EPCR thresholds ([official evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/evaluation.py)).

The official rules say participant data must be deleted within 30 days after the hackathon closes, submissions are released under CC BY, and publications must not add re-identifying information ([official rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/37e25dceda63ecec7c5b2ebeffd1ea0052ad886e/tabs/rules.py)).

**Recommendation.** Maintain two representations:

- a normalized internal representation for matching callers, annotation, and benchmarking; and
- an evaluator-facing representation produced by a deterministic final adapter and tested against the public evaluator.

Never allow a normalization step, contig-prefix conversion, or multiallelic decomposition to silently alter the final submitted allele tuple. Store the transformation history for every proposed allele.

### What the workflow must optimize

For each candidate pair, retain separate, inspectable components:

- gene–phenotype evidence;
- disease–phenotype evidence;
- inheritance compatibility;
- allele rarity and ancestry-aware frequency caveats;
- molecular consequence and predictor evidence;
- call confidence and read evidence for each allele;
- whether the two alleles are observed in cis, observed in trans, or unresolved;
- caller/aligner support provenance; and
- pair-level uncertainty converted to a calibrated EPCR.

Do not collapse these into an opaque model without retaining the inputs and component scores.

## 2. Small variants: aligners and callers that add real information

### Alignment choices

| Method | Architecture/status | Complementarity | Competition use |
|---|---|---|---|
| BWA-MEM2 2.3 | CPU-optimized BWA-MEM reimplementation. Its official repository says it should produce output identical to BWA-MEM 0.7.17 ([repository](https://github.com/bwa-mem2/bwa-mem2)). | **Not an independent aligner.** Faster BWA-MEM is operational diversity, not evidential diversity. | Sensible default linear alignment if already validated. Do not count BWA-MEM plus BWA-MEM2 as two votes. |
| DRAGMAP / open DRAGEN-GATK | Hash-table linear aligner derived from DRAGEN. Broad reports it as the aligner that replaces BWA-MEM in the open DRAGEN-GATK pipeline and emphasizes use of its tested masked hg38 reference ([official release](https://gatk.broadinstitute.org/hc/en-us/articles/4411716682011-Full-release-of-open-source-DRAGEN-GATK-1-0)). | A genuinely different mapper, though downstream GATK calling may overlap the primary caller. | Candidate-region or benchmark rescue; a full realignment is justified only if the primary mapping is suspicious. |
| vg Giraffe | Haplotype-aware pangenome graph mapper, not a linear-reference speed rewrite. The official documentation calls it a fast haplotype-aware short-read mapper; vg 1.76.0 was documented in July 2026 ([manual](https://github.com/vgteam/vg/wiki/vg-manpage)). | **Independent mapping architecture** and potentially helpful in polymorphic/divergent loci. | Late, targeted rescue or public-benchmark experiment. Current vg guidance says the ecosystem is more complex and less mature and recommends recent releases ([best practices](https://github.com/vgteam/vg/wiki/Giraffe-best-practices)); it is too operationally expensive for an unvalidated whole-genome default. |

**Recommendation.** Use one validated linear alignment as the backbone. Rerun a second aligner only for (a) public benchmark comparison, or (b) a small candidate interval where mapping quality, allele balance, soft clips, local repeats, pseudogene homology, or caller disagreement provides a reason. Alignment diversity without a gate multiplies compute and reconciliation work.

### Caller architectures and redundancy

| Caller | What is independent | 2026 status and terms | Recommended role |
|---|---|---|---|
| DeepVariant 1.10.0 | Converts local read evidence into pileup-like tensors and uses a learned neural classifier, unlike the explicit Bayesian assembly models below ([original method](https://www.nature.com/articles/nbt.4235)). | v1.10.0 released 2026-03-05; BSD-3-Clause code. The release mainly adds long-read phasing and infrastructure; use the explicit short-read `WGS` model ([official release](https://github.com/google/deepvariant/releases/tag/v1.10.0), [quick start](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-quick-start.md)). | Primary SNV/indel callset after GIAB calibration on matched read length, coverage, and alignment. Use CPU by default on this host; see GPU note below. |
| GATK HaplotypeCaller 4.7.0.0 | Detects active regions, performs local de novo assembly, aligns candidate haplotypes, calculates PairHMM likelihoods, and applies Bayesian genotyping ([official algorithm documentation](https://gatk.broadinstitute.org/hc/en-us/articles/4405451272731-HaplotypeCaller)). | GATK4 is Apache-2.0; 4.7.0.0 was the current release on 2026-08-25 ([official repository](https://github.com/broadinstitute/gatk), [releases](https://github.com/broadinstitute/gatk/releases)). | Independent rescue and concordance evidence, especially for complex local haplotypes. Retain allele likelihoods and assembly-region evidence. |
| Sentieon DNAseq | Sentieon describes DNAseq as a mathematically identical, optimized implementation of GATK algorithms ([official product page](https://www.sentieon.com/products/)). | Proprietary/commercial. | **Redundant with GATK HaplotypeCaller** for evidential purposes. It may substitute for speed, but do not ensemble DNAseq and GATK as independent callers. Sentieon DNAscope is a different proprietary ML caller, but adds licensing and transparency costs. |
| Strelka2 2.9.10 | Uses an indel error model, tiered haplotype modeling, and empirical variant rescoring distinct from DeepVariant and HaplotypeCaller ([primary paper](https://www.nature.com/articles/s41592-018-0051-x)). | Last tagged release 2018; the owner archived the repository on **2026-04-20**. Current `LICENSE.txt` is PolyForm Strict although legacy `COPYRIGHT.txt` contains GPL language; treat current use/redistribution as license-review required ([official archived repository](https://github.com/Illumina/strelka), [current license](https://github.com/Illumina/strelka/blob/v2.9.x/LICENSE.txt)). | Optional benchmark challenger, not a production pillar. No maintainer-endorsed active successor/fork was identified in the reviewed official material. DeepVariant/HaplotypeCaller are preferable maintained small-variant paths. |
| Octopus 0.7.4 | Unified haplotype-based Bayesian model supporting germline, somatic, arbitrary ploidy, and complex alleles ([primary paper](https://www.nature.com/articles/s41587-021-00861-3)). | MIT; last official release 0.7.4 in 2021, though repository activity continued ([official repository](https://github.com/luntergroup/octopus)). | Research-only third opinion on a bounded candidate interval or benchmark. Do not make an old release a critical workflow dependency without local validation. |

**Established:** DeepVariant, HaplotypeCaller, Strelka2, and Octopus are not merely different wrappers; their inference machinery differs. **Inference:** they still share the same reads, reference, and many mapping artifacts, so caller concordance is not independent biological replication.

### DeepVariant GPU decision for this host

DeepVariant’s official v1.10 quick start says only `call_variants` uses one GPU; `make_examples` and `postprocess_variants` remain CPU tasks ([official quick start](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-quick-start.md)). More importantly, the official FAQ reports that its Keras models occupy **16 GB of GPU memory** ([official FAQ](https://github.com/google/deepvariant/blob/r1.10/docs/FAQ.md#how-much-gpu-memory-is-needed-for-the-keras-models)). Therefore an RTX 5070 with 8 GB VRAM is **not a supported-safe assumption** for the v1.10 GPU container. Keep the CPU image as the reproducible default. A GPU attempt is acceptable only as a bounded public-test-data smoke benchmark that proves container/driver compatibility, peak VRAM, completion, output identity, and fallback behavior before any authorized case run; an out-of-memory result should end the experiment rather than trigger untracked model/runtime changes.

### Proposed small-variant consensus policy

Do **not** use majority vote. Use a union with provenance followed by calibrated evidence:

1. normalize a copy of each callset with a fixed reference and decompose only for internal matching;
2. reconcile equivalent haplotypes rather than exact VCF lines;
3. retain caller-specific genotype likelihood, depth, allele balance, strand/orientation, mapping quality, local complexity, and filter status;
4. classify support as primary-only, rescue-only, or concordant;
5. keep strong primary-only calls if read evidence is coherent; and
6. manually inspect both alleles of every top compound-het pair in a read viewer before submission.

Caller disagreement is a triage feature, not an automatic veto. Systematic mapping errors can fool every caller; an independent remapping is more useful than a fourth caller on the same BAM.

## 3. Compound heterozygosity and singleton phasing

### What short reads can and cannot establish

**Established.** WhatsHap performs read-based phasing. Its documentation says a read must cover at least two heterozygous calls and that even paired-end reads are only “somewhat helpful”; long reads work best ([official guide](https://github.com/whatshap/whatshap/blob/main/doc/guide.rst)). WhatsHap is MIT-licensed; its last numbered changelog release is 2.3 (2024), while development has continued ([repository/changelog](https://github.com/whatshap/whatshap/blob/main/CHANGES.rst)). HapCUT2 provides another read-backed maximum-likelihood haplotype assembly formulation ([primary paper](https://genome.cshlp.org/content/27/5/801), [official repository](https://github.com/vibansal/HapCUT2)).

**Inference.** With typical paired-end short reads, two heterozygous variants separated by more than the physical fragment span will have no molecule linking them. Without parents, long reads, linked reads, or informative relatives, their cis/trans relationship is normally unresolved. A statistical phaser can assign a phase, but a private pathogenic allele is precisely where population-based phase is least trustworthy.

SHAPEIT5 separates common-variant scaffold phasing from rare-variant phasing and reported accurate rare-variant phasing at biobank scale ([primary paper](https://www.nature.com/articles/s41588-023-01415-w), [official documentation](https://odelaneau.github.io/shapeit/)). However, its model is designed for large genotype datasets/reference panels, not proof of trans in one singleton; additionally, the official GitHub repository was disabled on the research date even though the documentation announced version 5.1.0. This is an operational and provenance blocker, not a foundation for a competition claim.

### Phase evidence vocabulary

For every candidate pair, report exactly one status:

- `trans_read_backed`: informative molecules support opposite haplotypes and no credible conflicting molecules;
- `cis_read_backed`: informative molecules support the same haplotype;
- `local_phase_conflict`: phase evidence is inconsistent or likely affected by mapping/duplicates;
- `unresolved_distance`: no molecule spans or bridges the alleles;
- `statistical_trans_support` / `statistical_cis_support`: population phasing only, with method, panel, ancestry fit, and posterior recorded; or
- `not_phaseable`: e.g. allele representation or locus quality makes the attempt invalid.

Only the first two are direct short-read molecular evidence. Statistical phase can reorder otherwise similar hypotheses but should not erase a strong pair or be described as proof.

### Pair construction without parents

**Recommendation.** Generate pair hypotheses after annotation rather than asking each caller to infer “compound het.” Within each gene and transcript-relevant disease model:

1. retain rare heterozygous alleles under an explicit maximum credible population frequency appropriate to a recessive disorder;
2. include coding, canonical splice, credible deep-splice, exon-level CNV/SV, and other allele classes in the same pair builder;
3. prevent a single complex haplotype from being double-counted as two causal alleles after decomposition;
4. reject or heavily penalize read-backed cis pairs;
5. preserve unresolved pairs, because lack of phase is expected in a singleton;
6. score allele quality separately, then combine at pair level; and
7. emit the exact two evaluator-facing alleles in one row.

Mixed-class pairs—such as SNV + exon deletion or coding allele + deep-intronic splice allele—are disproportionately important in a negative coding-only analysis and are a better use of extra compute than adding redundant SNV callers.

## 4. Orthogonal variant-class lanes

### SV and CNV

| Lane | Established capability/status | Singleton recommendation |
|---|---|---|
| DELLY 2.6.0 | Integrates paired-end, split-read, and read-depth evidence for deletions, insertions, duplications, inversions, translocations, and CNVs ([official repository](https://github.com/dellytools/delly)). BSD-3-Clause and actively released in August 2026. | First-pass breakpoint SV caller. Filter with official exclusion regions and inspect split/discordant reads. |
| CNVpytor 1.3.2 | Combines read depth with SNP/BAF information and supports local BAM/CRAM input ([primary paper](https://pubmed.ncbi.nlm.nih.gov/34817058/), [official repository](https://github.com/abyzovlab/CNVpytor)). MIT. | First-pass dosage/aneuploidy/allelic-imbalance lane. Use multiple bin sizes and require coherent segmentation; singleton calls lack cohort normalization. |
| Manta 1.6.0 | Paired/split-read plus local assembly, originally optimized for germline small sample sets ([official repository](https://github.com/Illumina/manta)). | Optional historical comparator only: the owner archived it on **2025-10-11**, and it is now under PolyForm Strict ([archive evidence](https://github.com/Illumina/manta/issues), [current license](https://github.com/Illumina/manta/blob/master/LICENSE.txt)). No official active successor/fork was identified in reviewed Illumina material; actively maintained DELLY2 is preferable for the first-pass singleton lane. |
| GATK-SV 1.1.x | Official single-sample workflow integrates multiple evidence types/callers against a precomputed reference panel ([single-sample documentation](https://gatk.broadinstitute.org/hc/en-us/articles/9022653744283-GATK-Best-Practices-for-Structural-Variation-Discovery-on-Single-Samples)); v1.1 adds single-sample filtering and a new panel ([release](https://github.com/broadinstitute/gatk-sv/releases/tag/v1.1)). Apache/BSD components, but individual caller terms still apply. | Comprehensive late-stage escalation if infrastructure and reference assets are ready. It is much heavier than DELLY + CNVpytor and primarily cloud/WDL oriented. |

**Recommendation.** Start with DELLY plus CNVpytor because their evidence is complementary (breakpoint versus dosage/BAF). Reconcile overlaps at the event/gene level; do not require identical breakpoints. Escalate to GATK-SV only after the first-pass lanes fail or identify an ambiguous high-value region. Explicitly look for a heterozygous exon deletion paired with an SNV/indel in the same gene.

### STR and repeat expansion

- **ExpansionHunter 5.0.0** is the current official release (2021-08-20). It uses a sequence-graph model to genotype catalogued repeats, including alleles longer than the read length ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6853681/)). Version 5 added multithreaded large-catalog support and reduced streaming memory use ([official v5.0.0 release](https://github.com/Illumina/ExpansionHunter/releases/tag/v5.0.0)). The repository remains available but its current official license is PolyForm Strict; code/package assets must be reviewed before use or redistribution ([official repository](https://github.com/Illumina/ExpansionHunter), [license](https://github.com/Illumina/ExpansionHunter/blob/master/LICENSE.txt)).
- **ExpansionHunter Denovo** detects repeat expansions without a predefined locus catalog and was designed for genome-wide discovery; its paper reports a catalog-free method and cohort/outlier workflows ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7187524/)). It is less informative in a singleton without matched controls and is now also PolyForm Strict ([repository](https://github.com/Illumina/ExpansionHunterDenovo)).
- **STRling** detects known and novel expansions from k-mer/soft-clipped evidence, is MIT-licensed, and had release 0.6.0 in 2025 ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380/), [official repository](https://github.com/quinlan-lab/STRling)).

**Recommendation.** Run a catalogued disease-repeat lane early if phenotype makes it plausible; use STRling as a license-friendly genome-wide screen. Treat EHdn singleton outliers as leads requiring locus visualization and independent confirmation, not as self-validating calls.

### Mobile-element insertions

**Established.** xTea combines clipped, discordant, and assembly/genotyping evidence across LINE1, Alu, SVA and other mobile elements ([primary paper](https://www.nature.com/articles/s41467-021-24041-8), [official repository](https://github.com/parklab/xTea)). The code is locally runnable, but its official license permits only internal non-commercial research at an academic/nonprofit institution and restricts distribution/derivative works ([official license](https://github.com/parklab/xTea/blob/master/LICENSE)).

**Recommendation.** Use xTea only after confirming the team/institution and hackathon use satisfy its license. It is a second-stage lane; inspect each top event because the official documentation warns that an orphan-transduction module has higher false-positive rates.

### Mitochondrial variation

**Established.** GATK’s mitochondrial short-variant workflow uses Mutect2, the standard and shifted mitochondrial references, merging, NUMT-oriented filters, and contamination/heteroplasmy filters ([official best practice](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels)).

**Recommendation.** Run the dedicated mtDNA workflow, not a diploid autosomal caller on chrM. Record heteroplasmy, strand/orientation evidence, depth, NUMT risk, and haplogroup context. This lane is biologically orthogonal but lower priority than a high-quality nuclear compound-het search given the evaluator’s public two-allele truth contract.

### Mosaicism and aneuploidy

**Established.** Mutect2 supports a single-sample mode and is designed to detect low-allele-fraction variants, but without a matched normal and panel of normals it cannot cleanly separate germline variants, artifacts, and true mosaicism ([official Mutect2 documentation](https://gatk.broadinstitute.org/hc/en-us/articles/27007991962907-Mutect2)). CNVpytor exposes chromosome-scale depth and BAF signals useful for copy-number/allelic imbalance.

**Inference.** The competition describes a real child but the public compound-heterozygous truth does not imply low-VAF somatic mosaicism. A whole-genome Mutect2 sweep is therefore a low-priority distraction unless allele balance, phenotype, or CNV/BAF evidence creates a specific mosaic hypothesis.

**Recommendation.** Perform a cheap chromosome-ploidy and BAF screen. Use Mutect2 only in gated candidate regions or as a late safety lane, with artifact flags and tissue limitation stated explicitly.

### Deep intronic, splice, and regulatory variation

**Established.** SpliceAI is a deep residual sequence model for splice gain/loss prediction ([primary paper](https://doi.org/10.1016/j.cell.2018.12.015)); its repository was archived in April 2026, code is PolyForm Strict, and bundled models are CC BY-NC 4.0 for academic/non-commercial use ([official license](https://github.com/Illumina/SpliceAI/blob/master/LICENSE)). Pangolin predicts tissue-aware splice-site strength using a separately trained deep model, runs locally on VCF/CSV, and is GPL-3.0 ([primary paper](https://doi.org/10.1186/s13059-022-02664-4), [official repository](https://github.com/tkzeng/Pangolin)).

Genomiser extends phenotype-aware prioritization into noncoding sequence using ReMM and optional CADD data; Exomiser 15.1.0 officially documents local ReMM/CADD configuration and recommends 12–16 GB heap for genome analysis ([official installation guide](https://exomiser.readthedocs.io/en/stable/installation.html)). The 100,000 Genomes pilot reported noncoding diagnoses requiring targeted whole-genome analysis and functional validation, including variants paired with loss-of-function alleles ([primary clinical study](https://www.nejm.org/doi/full/10.1056/NEJMoa2035790)).

**Recommendation.** Avoid an unfiltered genome-wide “regulatory AI score” hunt. Prioritize:

1. canonical splice and splice-region alleles;
2. deep intronic SpliceAI/Pangolin candidates in top phenotype genes;
3. a deep intronic candidate as the second allele to a strong coding/SV allele; and only then
4. promoter/UTR/enhancer candidates with credible gene linkage.

SpliceAI and Pangolin scores are supporting predictions, not pathogenicity proof. Preserve distance to predicted cryptic splice site, affected transcript, tissue expression, transcript choice, and cross-model agreement. Functional RNA/minigene confirmation is outside the competition workflow but should be named as the downstream validation need.

## 5. Phenotype, gene, disease, and variant prioritization

### Primary local tools

| Tool | Model and transparency | Current status / privacy / terms | Best use here |
|---|---|---|---|
| Exomiser / Genomiser 15.1.0 + 2602 data | Combines variant frequency/quality/pathogenicity, inheritance models, human disease and model-organism phenotype similarity; emits gene/variant rankings and decomposable results ([official overview](https://github.com/exomiser/Exomiser), [method](https://pmc.ncbi.nlm.nih.gov/articles/PMC5467691/)). | **Application:** v15.1.0 is the current tagged release (2026-06-09), Java 21, AGPL-3.0 ([release](https://github.com/exomiser/Exomiser/releases/tag/15.1.0)). **Data:** 2602 is the newest official data announcement (2026-05-13), with gnomAD v4.1, dbNSFP 4.5a, ClinVar 2026-02-08 and updated phenotype sources ([official 2602 announcement](https://github.com/exomiser/Exomiser/discussions/640)). The announcement predates 15.1.0 and explicitly lists compatibility through 15.0.0/14.x; development docs show 15.1.1 + 2602 but do not establish a release. | Main automated phenotype + inheritance-aware prioritizer. Pin 15.1.0 + 2602 only after the official test analysis and a public fixture pass; record this locally validated combination rather than claiming it from the release page. Fully offline after download. |
| LIRICAL 2.4.1 | Calculates interpretable likelihood ratios and shows the contribution of each HPO feature and genotype to each disease ([primary paper](https://doi.org/10.1016/j.ajhg.2020.06.021), [algorithm docs](https://thejacksonlaboratory.github.io/LIRICAL/stable/explanations.html)). | v2.4.1 released 2026-06-18; local. License is non-commercial, non-transferable/non-sublicensable, and prohibits modification/derivative works ([official license](https://github.com/TheJacksonLaboratory/LIRICAL/blob/master/LICENSE)). | Independent disease-level phenotype audit and explanation. Do not redistribute it inside a public workflow image without legal review. |
| Phen2Gene 1.2.3 | Phenotype-only gene ranking from HPO; no variant or inheritance model ([primary paper](https://doi.org/10.1093/nargab/lqaa032)). | MIT/local, but latest release is from 2021 ([official repository](https://github.com/WGLab/Phen2Gene)). | Lightweight independent phenotype-only sensitivity check, not the primary ranker. |

### Rank fusion without opacity

**Recommendation.** Do not average tool ranks. Build a transparent pair evidence table and use a small, prespecified fusion rule. For example:

`pair_score = phenotype_component + inheritance_component + allele1_component + allele2_component + phase_component + orthogonal_evidence_component - artifact_penalties`

Each component should have a documented scale, missingness behavior, and cap. Preserve both Exomiser and LIRICAL outputs as features, because one is a joint gene/variant prioritizer while the other is a disease-level likelihood-ratio explanation. A phenotype-only method can expose rank sensitivity to the variant filters.

Run phenotype sensitivity analyses:

- all provided HPO terms;
- high-specificity terms only;
- one-term-at-a-time deletion;
- with and without negated terms where supported;
- parent/ancestor term collapse to test ontology redundancy; and
- phenotype-only versus phenotype-plus-variant modes.

A pair that disappears after removing one uncertain HPO term is less stable than one supported across perturbations. Report rank ranges, not just the best rank.

### Privacy and offline operation

**Recommendation.** Keep the proband VCF, HPO profile, BAM/CRAM, candidate list, and generated explanations local. Download public ontologies/models/resources once, checksum them, and run offline. Do not send row-level data to hosted LLMs, web prioritizers, notebooks, analytics endpoints, or telemetry-enabled services without explicit organizer approval and a documented data-processing agreement. The public challenge rules impose privacy and deletion duties; a public web interface saying it deletes uploads is not equivalent to authorization.

## 6. Annotation and pathogenicity resource snapshot

| Resource | Established use | Local availability and 2026 snapshot | License/operational caveat |
|---|---|---|---|
| Ensembl VEP | Transcript consequence engine and plugin framework. | VEP 116.1 released 2026-08-04; caches and FASTA support fully offline operation ([official repository/releases](https://github.com/Ensembl/ensembl-vep/releases), [download docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)). Apache-2.0 code. | Cache, FASTA, transcript set, MANE version, and every plugin data file must be pinned together. Plugin resources retain their own terms. |
| gnomAD | Population allele frequency, coverage/AN, constraint, CNV/SV and repeat context. | gnomAD v4.1.1 released 2026-03-30 with updated constraint/LOFTEE flags and VEP 115 annotations ([official release](https://gnomad.broadinstitute.org/news/2026-03-gnomad-v4-1-1/)). Downloadable Hail/VCF resources are local; v4 is GRCh38. | Do not equate absence with zero frequency when AN/coverage is low. Preserve ancestry groups, FAF, filters, exome/genome discordance, and version. Check gnomAD terms before redistribution. |
| ClinVar | Submitted variant–condition assertions, review status, conflicts, evidence links. | Web/API/FTP updated weekly; monthly releases are archived ([official release cycle](https://www.ncbi.nlm.nih.gov/clinvar/docs/release_cycle/), [download guide](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)). | VCF is partial (<10 kb simple mapped alleles); use variant-centric XML/TSV for full assertion context. Never collapse conflicts or star/review status into a binary “ClinVar pathogenic” flag. |
| ClinGen | Expert gene–disease validity, dosage sensitivity, actionability, and expert curation. | Curated downloads/APIs can be mirrored locally. ClinGen releases curated content under CC0 and requests attribution/date accessed ([official terms](https://www.clinicalgenome.org/docs/terms-of-use/)). | Use gene–disease validity and dosage evidence separately from variant assertions. External linked resources may have different terms. |
| LOFTEE | VEP plugin classifying putative LoF alleles as high/low confidence with transcript and rescue filters. | Local VEP plugin; MIT. GRCh38 requires the GRCh38 branch, not master ([official README](https://github.com/konradjk/loftee/blob/master/README.md)). gnomAD v4.1.1 updated LOFTEE behavior/flags. | Pin branch, VEP/transcript version, ancestral FASTA, and conservation database. “HC LoF” is a consequence-quality filter, not proof of disease causality. |
| SpliceAI | Long-context neural splice gain/loss prediction. | Local code/models technically available; repository archived 2026. | PolyForm Strict code; CC BY-NC 4.0 models; commercial use needs a license ([official license](https://github.com/Illumina/SpliceAI/blob/master/LICENSE)). Bundled annotations are old GENCODE canonical sets unless replaced deliberately. |
| Pangolin | Tissue-aware splice-site strength prediction with an independent deep model. | Local CLI; latest tagged model/software 1.0.1 (2022); GPL-3.0 ([repository](https://github.com/tkzeng/Pangolin)). | Requires a transcript database; simple-variant focus. Use as corroboration/sensitivity evidence. |
| CADD | Genome-wide deleteriousness score integrating conservation and functional annotations. | CADD score model v1.7; script release 1.7.3 (2025). Offline scoring needs roughly 100 GB–1 TB disk and at least 12 GB RAM; pre-scored files are simpler ([official repository](https://github.com/kircherlab/CADD-scripts), [v1.7 notes](https://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.7.pdf)). | Official license grants broad use to non-commercial users/licensees; commercial use requires a license. It is a general rank feature, not calibrated Mendelian probability. |
| REVEL | Random-forest ensemble for rare missense variants, combining multiple functional/conservation predictors ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/)). | Precomputed v1.3 (May 2021) includes transcript IDs and GRCh37/38 positions; local VEP plugin supported ([official VEP plugin docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)). | Ensembl documents REVEL as non-commercial. It is missense-only and includes constituent predictors also present in other ensembles; correlated scores are not independent evidence. |
| AlphaMissense | Proteome-wide missense effect classification from protein sequence/structure-informed learning ([primary paper](https://www.science.org/doi/10.1126/science.adg7492)). | Precomputed human predictions can be downloaded and used locally. Reference code repository was archived; weights are not released. Code Apache-2.0, predictions CC BY 4.0 ([official repository](https://github.com/google-deepmind/alphamissense)). | Official disclaimer says it is not validated/approved for clinical use. Transcript/isoform mapping and missing predictions must be explicit. Do not treat its class as an ACMG verdict. |

### Annotation rules that matter more than another score

1. Pin GRCh38 FASTA digest, contig dictionary, VEP/cache release, GENCODE/Ensembl and MANE releases.
2. Retain all relevant transcripts, then identify MANE Select/Plus Clinical and disease-relevant alternatives; do not choose the most severe consequence across unrelated transcripts without showing which transcript produced it.
3. Keep raw allele frequency counts and AN, not only AF.
4. Treat predictors from the same training labels/features as correlated. CADD, REVEL, AlphaMissense, SpliceAI, and Pangolin answer different questions; none proves pathogenicity.
5. Pin ClinVar monthly archive for reproducibility while optionally checking the current weekly record during final manual review.
6. Record license, source URL, version/date, checksum, assembly, contig convention, and redistribution status for every resource.

## 7. Benchmarking and calibration

### Public truth layers

**Small variants and the v5.0q reference boundary.** NIST now lists HG002 v5.0q as the current HG002 small-variant/SV benchmark and deprecates v4.2.1 for HG002. However, NIST also states that v5.0q is constructed from the v1.1 **T2T-HG002 personal assembly**, while distributing projected representations for GRCh37, GRCh38, and T2T-CHM13 ([current NIST GIAB page](https://www.nist.gov/programs-projects/genome-bottle), [NIST FAQ](https://www.nist.gov/programs-projects/faqs-genome-bottle)). It is therefore not a direct like-for-like truth construction for this conventional GRCh38 short-read workflow, even when using its GRCh38 projection.

**Recommendation:** retain GIAB HG002 **v4.2.1 on the exact GRCh38 reference bundle** as the primary reference-matched calibration, while labeling it deprecated by NIST. It covers substantially more segmental duplication and low-mappability sequence than v3.3.2 but still excludes difficult sequence and SVs ([primary v4.2.1 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9706577/)). Use current GIAB **v3.6 stratifications** for GRCh38 context analysis ([official NIST page](https://www.nist.gov/programs-projects/genome-bottle), [stratification repository](https://github.com/usnistgov/giab-stratifications)). Add v5.0q’s GRCh38 projection as a separate sensitivity benchmark only after following its README and verifying FASTA identity, projection assumptions, confident regions, and comparison-tool support. Do not mix v4.2.1 and v5.0q counts into one headline metric.

**Challenging medically relevant genes.** Run the separate GIAB CMRG v1.00 small-variant and SV benchmarks across 273 difficult medically relevant autosomal genes. The authors explicitly recommend using CMRG alongside, not merged into, v4.2.1 because the benchmarks have different scope and construction ([primary paper and data links](https://www.nature.com/articles/s41587-021-01158-1), [NIST catalog](https://catalog.data.gov/dataset/challenging-medically-relevant-genes-benchmark-set)).

**SV.** Use the GIAB HG002 v0.6 Tier 1 SV benchmark and CMRG SV truth, reporting their restricted regions separately. Do not advertise genome-wide SV sensitivity from a Tier 1 subset.

**Mosaic.** If a mosaic lane is retained, use NIST’s public HG002 mosaic benchmark as a limited positive set, while stating its chemistry/depth and truth-design limitations ([official NIST repository](https://github.com/usnistgov/giab-HG002-mosaic-benchmark)).

### Comparison and normalization

**Established.** GA4GH best practices require truth VCF, query VCF, confident regions, haplotype-aware comparison, and stratification. Simple left-shifting/trimming cannot reconcile all equivalent complex representations; hap.py with vcfeval was the reference approach ([GA4GH primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/), [official tools repository](https://github.com/ga4gh/benchmarking-tools)). vcfdist adds local-phase-aware distance matching and phase-error evaluation and is useful as a sensitivity comparator ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10710436/), [official repository](https://github.com/TimD1/vcfdist)).

**Recommendation.** For each caller/alignment configuration:

- report precision, recall, F1 and genotype errors for SNVs and indels separately;
- stratify by low complexity, homopolymers, segmental duplications, mappability, GC, coding regions, difficult medically relevant genes, and GIAB compound-het/complex regions using official stratifications ([GIAB stratifications](https://github.com/usnistgov/giab-stratifications));
- report outside-confident-region call counts separately, never as presumed errors or successes;
- preserve the pre-normalized callset and a transformation manifest; and
- test that the final evaluator adapter recreates the literal expected tuple after internal haplotype reconciliation.

### Synthetic compound-het tests

Public truth sets test calling, but not the full pair-ranking problem. Add two synthetic layers:

1. **VCF-level ranking fixtures:** create public, non-patient mini-cases with two truth alleles in the same gene and distractors. Vary allele consequences, distances, cis/trans metadata, missing phase, frequency, transcript, one allele being an SV, and HPO specificity. These are deterministic unit/e2e tests of pairing, rank fusion, EPCR, and submission serialization—not caller benchmarks.
2. **Read-level spike-ins:** use a public GIAB BAM and insert paired alleles on controlled haplotypes for a small design matrix. BAMSurgeon can add SNVs/indels/SVs and its own manual warns that haplotype-sensitive experiments should phase first and spike per haplotype ([official repository/manual](https://github.com/adamewing/bamsurgeon)). VarSim can construct diploid synthetic genomes and reads ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4410653/)).

**Caveat.** Simulation does not reproduce every mapping, library, ancestry, repeat, or sequencing artifact. The BAMSurgeon maintainers explicitly recommend validated truth data when suitable. Therefore use GIAB for caller performance and synthetic pairs for controlled failure discovery, not as evidence of clinical accuracy.

### Ranking and EPCR calibration

Measure the full pipeline on public/synthetic cases with:

- pair Recall@1, @3, @5, and @10;
- fraction with both causal alleles retained before ranking;
- allele-level recall versus pair-level recall;
- median reciprocal rank of the true pair;
- top-k candidate count by gene and by unique allele;
- calibration curve, Brier score, and log loss for any EPCR-like probability;
- F-max under the **official evaluator implementation**; and
- rank stability under tool/resource/phenotype perturbations.

Never call a heuristic score a probability. Fit a monotonic calibration (for example isotonic regression) only if there are enough held-out public/synthetic cases; otherwise use conservative EPCR bands with documented meaning. Avoid tied EPCRs because the evaluator breaks ties by file order.

### Sensitivity analysis matrix

At minimum, compare:

- DeepVariant alone versus DeepVariant union HaplotypeCaller;
- default alignment versus gated alternative alignment in difficult strata;
- ClinVar/gnomAD/annotation release snapshots;
- phenotype term subsets and leave-one-out profiles;
- phase unknown versus statistical phase versus read-backed phase;
- frequency thresholds and quality thresholds;
- coding-only versus coding + splice-targeted versus all noncoding;
- with/without predictor families (missense and splice separately); and
- candidate pair ranks before and after SV/CNV/repeat lanes.

Promote a workflow change only if it improves held-out pair recall/rank stability without unacceptable false-pair inflation. A one-off improvement in one simulated case is not sufficient.

## 8. Staged, compute-efficient execution with review checkpoints

### Stage 0 — contract and provenance (cheap, mandatory)

- fingerprint FASTA, BAM/CRAM, index, read groups, sample ID, contigs, and all public resources;
- record tool/container/Pixi-lock digests and commands;
- check BAM sort/index, duplicates, coverage distribution, contamination/sex/ploidy indicators, insert size, and chrM coverage;
- prove the submission adapter against the immutable public evaluator with synthetic rows.

**Inspect:** resolve reference, contig, or sample-identity inconsistencies before
interpreting downstream results. Do not “fix” mismatches silently.

### Stage 1 — primary coding/splice SNV/indel pass (highest yield)

- DeepVariant WGS on the validated alignment;
- fixed internal normalization and VEP/LOFTEE annotation;
- local gnomAD, ClinVar/ClinGen, missense and splice annotation;
- Exomiser recessive/all-MOI analysis plus LIRICAL phenotype audit;
- gene-level pair generation including canonical splice and near-splice alleles.

**Inspect:** review the top pair table, not only top variants. Confirm each pair
has two distinct alleles, consistent zygosity, a valid transcript/gene
relationship, a plausible recessive disease model, and no representation
duplication.

### Stage 2 — independent small-variant rescue

- HaplotypeCaller on the whole genome if compute is available, otherwise all top phenotype genes plus flanks and difficult loci;
- haplotype-aware reconciliation rather than line intersection;
- WhatsHap/HapCUT2 read-backed phasing for candidate pairs;
- read-view inspection of both alleles and local haplotype complexity.

**Inspect:** use DRAGMAP or Giraffe where mapping evidence is suspect or in a
bounded public benchmark. Treat a rescue-only candidate as provisional until
its reads and phenotype/inheritance evidence are coherent.

### Stage 3 — orthogonal variant classes

- DELLY + CNVpytor;
- disease-repeat catalog plus STRling;
- GATK mitochondrial workflow;
- xTea only if license and phenotype/evidence gate allow;
- cross-class pair building (SNV + CNV/SV/splice/MEI).

**Inspect:** manually review event evidence and gene overlap. Do not translate
imprecise SVs into exact small-variant coordinates.

### Stage 4 — targeted noncoding rescue

- SpliceAI and Pangolin in phenotype-supported genes and around strong first alleles;
- Genomiser/ReMM as a broader sensitivity lane;
- promoter/UTR/regulatory review only for top genes with a credible gene-regulatory link.

**Inspect:** record transcript-aware splice/regulatory evidence, rarity,
coherent reads, and whether a plausible second allele exists. Keep functional
validation explicitly unresolved where it has not been performed.

### Stage 5 — final pair adjudication and manual submission

- rerun all top-pair coordinates directly against reference and alignment;
- verify evaluator-facing contig, position, REF, ALT, and paired row;
- deduplicate alleles and pairs;
- assign calibrated, strictly ordered EPCRs;
- run the public evaluator locally on format-only/synthetic checks;
- produce an evidence packet for every row and a deletion manifest for challenge close;
- prepare and upload the final file manually, outside the analysis workflow.

During exploration, rerun and replace derived artifacts freely. For results that
will be published, retain the exact Pixi lock, configuration, command log,
source versions, and evidence packet used for the reported result.

## 9. What is high yield for this competition

### Highest yield

1. **Exact pair construction and serialization.** Full credit requires the exact pair in one row.
2. **Phenotype/inheritance-aware gene ranking.** Exomiser plus a LIRICAL audit directly targets rare-disease ranking while remaining explainable.
3. **One strong caller plus one independent rescue caller.** DeepVariant + HaplotypeCaller captures meaningful architectural diversity.
4. **Candidate-level read review and phase status.** It prevents clean-looking but cis, duplicated, or mapping-artifact pairs from dominating.
5. **Mixed-class second alleles.** CNV/SV or deep-splice alleles paired with a strong coding allele are a biologically plausible failure mode of coding-only search.
6. **Public benchmark and synthetic pair e2e tests.** They find broken normalization, pairing, ranking, and evaluator transformations before submission.
7. **Rank/EPCR stability.** The evaluator rewards ordering and threshold behavior; calibration is part of the model, not presentation polish.

### Medium yield, gated

- DELLY + CNVpytor and repeat screens;
- dedicated mitochondrial calling;
- targeted alternative alignment in difficult loci;
- Pangolin/SpliceAI in top genes;
- GATK-SV comprehensive single-sample escalation;
- mobile-element calling when phenotype or local evidence supports it.

### Low value or actively risky

- BWA-MEM and BWA-MEM2 counted as independent evidence;
- GATK HaplotypeCaller and Sentieon DNAseq counted as independent evidence;
- four or five callers combined by majority vote;
- whole-genome graph remapping before the linear workflow is benchmarked;
- unrestricted whole-genome regulatory scoring without phenotype/gene gates;
- statistical phase described as proof of trans;
- treating multiple correlated missense ensemble scores as separate votes;
- uploading case data to hosted LLMs/web tools without explicit authorization;
- adopting archived/restrictively licensed Illumina tools without a license audit;
- optimizing only allele-level recall when the evaluator scores a pair hypothesis;
- tuning thresholds on hidden leaderboard feedback or any private answer-derived signal.

## 10. Minimum transparent evidence record

For every submitted pair, the internal record should include:

- immutable pair ID and exact evaluator tuple;
- normalized internal alleles and transformation trace;
- genome/reference/contig/transcript identifiers;
- each caller and aligner source with version, filter, likelihood/quality, depth, allele balance, mapping and strand/orientation evidence;
- gene, disease, and mode-of-inheritance rationale;
- HPO inputs and per-tool phenotype contributions;
- gnomAD AC/AN/AF/FAF by group, coverage/quality flags;
- ClinVar assertions, review status, conflicts, condition match, and access release;
- ClinGen gene–disease/dosage evidence;
- consequence, LOFTEE, missense and splice scores with transcript/model versions;
- phase status and raw supporting/conflicting molecule counts;
- SV/CNV/repeat/mt/MEI evidence where relevant;
- manual review decision and rationale;
- rank under every sensitivity run, not only final rank;
- EPCR derivation/calibration bin; and
- final file checksum and evaluator dry-run result.

This is enough to reconstruct why a pair was promoted without exposing it outside the authorized local environment.

## 11. Unresolved blockers and decisions required before implementation

1. **Competition data-processing boundary:** the public rules establish privacy/deletion duties but do not, in the reviewed source, authorize particular hosted AI or external analysis services. Keep all patient-level processing offline unless organizers explicitly clarify it.
2. **Reference bundle:** the exact FASTA build, decoys/alt contigs, contig naming, and original alignment recipe must be established from authorized metadata before choosing DRAGMAP/Giraffe or interpreting literal evaluator coordinates.
3. **License fit:** SpliceAI, Strelka2, Manta, ExpansionHunter/EHdn, LIRICAL, xTea, CADD, and REVEL have noncommercial/restrictive or mixed terms. Confirm that team status, prize competition use, public code/report distribution, containers, and result redistribution are permitted. Do not redistribute restricted code/models/data inside the public repository.
4. **Resource redistribution:** Exomiser application code is AGPL-3.0, but its downloadable data aggregate multiple upstream resources. Create a resource-by-resource bill of materials before publishing an environment or cache.
5. **Statistical phasing operations:** SHAPEIT5’s documentation is live but its GitHub repository was disabled as of the research date. Do not depend on it until source provenance and a reproducible version are restored; in any case, singleton rare-variant statistical phase remains supporting evidence only.
6. **Mosaic lane priority:** absent a specific low-VAF or chromosomal-imbalance signal, a full mosaic caller sweep is not justified by the public two-allele competition objective.

## Bottom line

The state-of-the-art workflow for this task is a transparent decision system, not a software collection. DeepVariant plus HaplotypeCaller, Exomiser plus LIRICAL, candidate read-backed phasing, DELLY plus CNVpytor, a repeat/mt safety net, and gene-focused splice rescue cover the important independent evidence classes. The decisive engineering is in allele equivalence, pair construction, cross-class pairing, rank calibration, and exact evaluator serialization. Those are also the parts that can be benchmarked honestly without touching or learning from the private case answer.
