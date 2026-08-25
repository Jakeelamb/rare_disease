# External reference resources

Large public reference databases are not committed. Each workflow run records
their release and checksums in its provenance manifest.

## Ensembl VEP

Install every locked Pixi environment, then use the isolated VEP environment
for the Ensembl 116 GRCh38 cache and FASTA:

```bash
pixi install --all --locked
pixi run -e vep vep_install -a cf -s homo_sapiens -y GRCh38 \
  -c resources/vep --CONVERT --USE_HTTPS_PROTO
```

The workflow invokes the VEP environment's Perl and launcher by explicit path,
so system Perl cannot silently take precedence.

For a resumable, checksum-gated cache-only install, the workflow's current
resource can also be recreated directly:

```bash
mkdir -p resources/vep/tmp
aria2c --continue=true --max-connection-per-server=16 --split=16 \
  --min-split-size=16M --file-allocation=none --dir=resources/vep/tmp \
  --out=homo_sapiens_vep_116_GRCh38.tar.gz \
  'https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz'
curl --fail --location \
  'https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/CHECKSUMS' \
  --output resources/vep/tmp/CHECKSUMS
sum resources/vep/tmp/homo_sapiens_vep_116_GRCh38.tar.gz
# Expected upstream sum fields: 56036 26996736
tar --extract --gzip \
  --file resources/vep/tmp/homo_sapiens_vep_116_GRCh38.tar.gz \
  --directory resources/vep
```

The agnostic coding lane bounds VEP work with the release-matched Ensembl 116
GTF. It retains all canonical protein-coding exons plus 20-bp splice flanks;
this is an efficiency boundary, not a whole-genome claim:

```bash
curl --fail --location \
  'https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz' \
  --output resources/vep/Homo_sapiens.GRCh38.116.gtf.gz
sum resources/vep/Homo_sapiens.GRCh38.116.gtf.gz
# Expected upstream sum fields: 49151 137815
gzip -dc resources/vep/Homo_sapiens.GRCh38.116.gtf.gz \
| awk 'BEGIN {OFS="\t"} \
  $1 ~ /^([0-9]+|X|Y|MT)$/ && $3 == "exon" && \
  $0 ~ /gene_biotype "protein_coding"/ \
  {chrom=($1 == "MT" ? "M" : $1); start=$4-21; if (start < 0) start=0; \
   print chrom,start,$5+20}' \
| sort -k1,1V -k2,2n -k3,3n \
| bedtools merge -i - \
> resources/vep/Ensembl.GRCh38.116.protein_coding_exons_plus20.bed
sha256sum resources/vep/Ensembl.GRCh38.116.protein_coding_exons_plus20.bed
# Expected SHA-256: 3110a828f7a57e3fdbad06c00906ce2fa89f63c8147ed88110aed2027531a307
```

## Human Phenotype Ontology

The phenotype lane uses the public HPO `v2026-06-23` release:

```bash
mkdir -p resources/hpo
curl --fail --location \
  'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-06-23/hp.obo' \
  --output resources/hpo/hp.obo
curl --fail --location \
  'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-06-23/genes_to_phenotype.txt' \
  --output resources/hpo/genes_to_phenotype.txt
sha256sum --check <<'CHECKSUMS'
a5092cbdf605f568403cf7380d9173014015692433b2cc631bc5c1b053876b1b  resources/hpo/hp.obo
26cb7ee00c73b5777f6e5ad43323c941e1fcef1d191592f332d7929f3ea1ab3f  resources/hpo/genes_to_phenotype.txt
CHECKSUMS
```

## Alignment reference

The source VCF names
`GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`.
A superficially similar Broad masked reference was tested and rejected: it has
3,366 contigs, including ALT loci, versus the source VCF's 2,580 contigs.
`reference_manifest.yaml` preserves that negative result.

Rebuild the named source contract from two primary NCBI resources: the GIAB v2
masked no-ALT primary reference plus NCBI's no-ALT + hs38d1 combined reference.
The Python builder verifies primary order and length before substitution:

```bash
mkdir -p resources/reference
curl --fail --location --continue-at - \
  'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz' \
  --output resources/reference/masked_primary_v2.fasta.gz
curl --fail --location --continue-at - \
  'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz' \
  --output resources/reference/no_alt_plus_hs38d1.fna.gz
printf '%s  %s\n' \
  a056c57649f3c9964c68aead3849bbf8 \
  resources/reference/no_alt_plus_hs38d1.fna.gz | md5sum --check
pixi run -e wgs python scripts/build_source_reference.py \
  --masked-primary resources/reference/masked_primary_v2.fasta.gz \
  --combined resources/reference/no_alt_plus_hs38d1.fna.gz \
  --output resources/reference/source_equivalent.no_chr.fasta
pixi run -e wgs samtools faidx resources/reference/source_equivalent.no_chr.fasta
pixi run -e wgs bcftools norm --check-ref e \
  --fasta-ref resources/reference/source_equivalent.no_chr.fasta \
  --output-type u data/WGS_EX2312012_HGWCNDSX7.vcf.gz >/dev/null
```

Only after the contig and REF checks pass may this FASTA be used for orthogonal
alignment. Do not substitute Ensembl's annotation FASTA: although primary
assembly coordinates agree, its contig/masking contract is not the one recorded
by the provided caller.
