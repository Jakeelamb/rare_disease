"""Opt-in WGS recall and an independent caller-comparison lane."""

wgs_root = f"{private_root}/wgs"
wgs_reference = config["references"]["alignment_fasta"]
wgs_lane_bams = expand(f"{wgs_root}/alignment/{{lane}}.name.bam", lane=lanes)
deepvariant_image = config["containers"]["deepvariant"]
read_audit = config["read_audit"]
phase_contigs = tuple([str(chromosome) for chromosome in range(1, 23)] + ["X", "Y", "M"])
phase_shard_vcfs = expand(
    f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.phased.vcf.gz",
    phase_contig=phase_contigs,
)
phase_shard_tbis = expand(
    f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.phased.vcf.gz.tbi",
    phase_contig=phase_contigs,
)


rule wgs_recall:
    input:
        f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz",
        f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz.tbi",
        f"{wgs_root}/deepvariant/PROBAND01.phased.vcf.gz",
        f"{wgs_root}/deepvariant/PROBAND01.phased.vcf.gz.tbi",
        f"{wgs_root}/qc/PROBAND01.whatshap.tsv",
        f"{wgs_root}/qc/PROBAND01.flagstat.txt",
        f"{wgs_root}/qc/leading_candidate_reads.json",
        f"{wgs_root}/qc/leading_candidate_reads.review.md",


rule wgs_comparison:
    input:
        f"{wgs_root}/qc/supplied_vs_deepvariant.bcftools_stats.txt",
        f"{wgs_root}/qc/agnostic_top30_reads.review.md",
        f"{wgs_root}/ranking/deepvariant_coding.review.md",
        f"{wgs_root}/ranking/caller_comparison.md",
        f"{wgs_root}/adjudication/leading_candidate.json",
        f"{wgs_root}/adjudication/leading_candidate.review.md",


rule index_wgs_reference:
    input:
        fasta=wgs_reference,
        fai=f"{wgs_reference}.fai",
    output:
        f"{wgs_reference}.0123",
        f"{wgs_reference}.amb",
        f"{wgs_reference}.ann",
        f"{wgs_reference}.bwt.2bit.64",
        f"{wgs_reference}.pac",
    threads: config["threads"]["alignment"]
    log:
        f"{private_root}/logs/bwa_index.log",
    params:
        bwa=str(wgs_bin / "bwa-mem2"),
    shell:
        "{params.bwa:q} index {input.fasta:q} > {log:q} 2>&1"


rule align_lane_name_sorted:
    input:
        r1=lambda wildcards: read_pairs[wildcards.lane][0],
        r2=lambda wildcards: read_pairs[wildcards.lane][1],
        reference=wgs_reference,
        indexes=rules.index_wgs_reference.output,
    output:
        temp(f"{wgs_root}/alignment/{{lane}}.name.bam"),
    threads: config["threads"]["alignment"]
    resources:
        mem_mb=60000,
    params:
        read_group=lambda wildcards: (
            f"@RG\\tID:{wildcards.lane}\\tSM:PROBAND01\\tLB:WGS\\tPL:ILLUMINA"
        ),
        bwa=str(wgs_bin / "bwa-mem2"),
        samtools=str(wgs_bin / "samtools"),
    log:
        f"{private_root}/logs/align_{{lane}}.log",
    wildcard_constraints:
        lane=r"L\d{3}",
    shell:
        """
        {params.bwa:q} mem -t {threads} -R {params.read_group:q} {input.reference:q} \
          {input.r1:q} {input.r2:q} 2> {log:q} \
        | {params.samtools:q} sort -n -@ 4 -O BAM -o {output:q} - 2>> {log:q}
        """


rule merge_name_sorted_lanes:
    input:
        wgs_lane_bams,
    output:
        temp(f"{wgs_root}/alignment/PROBAND01.merged.name.bam"),
    threads: config["threads"]["alignment"]
    log:
        f"{private_root}/logs/merge_name_sorted.log",
    params:
        samtools=str(wgs_bin / "samtools"),
    shell:
        """
        {params.samtools:q} merge -n -@ {threads} -O BAM -o {output:q} \
          {input:q} 2> {log:q}
        """


rule fixmate_coordinate_sort:
    input:
        f"{wgs_root}/alignment/PROBAND01.merged.name.bam",
    output:
        temp(f"{wgs_root}/alignment/PROBAND01.fixmate.coord.bam"),
    threads: config["threads"]["alignment"]
    resources:
        mem_mb=60000,
    log:
        f"{private_root}/logs/fixmate_sort.log",
    params:
        samtools=str(wgs_bin / "samtools"),
    shell:
        """
        {params.samtools:q} fixmate -m -@ {threads} {input:q} - 2> {log:q} \
        | {params.samtools:q} sort -@ {threads} -O BAM -o {output:q} - 2>> {log:q}
        """


rule mark_duplicates:
    input:
        f"{wgs_root}/alignment/PROBAND01.fixmate.coord.bam",
    output:
        bam=f"{wgs_root}/alignment/PROBAND01.markdup.bam",
        bai=f"{wgs_root}/alignment/PROBAND01.markdup.bam.bai",
        flagstat=f"{wgs_root}/qc/PROBAND01.flagstat.txt",
    threads: config["threads"]["alignment"]
    log:
        f"{private_root}/logs/markdup.log",
    params:
        samtools=str(wgs_bin / "samtools"),
    shell:
        """
        {params.samtools:q} markdup -s -@ {threads} {input:q} {output.bam:q} 2> {log:q}
        {params.samtools:q} index -@ {threads} {output.bam:q} {output.bai:q} 2>> {log:q}
        {params.samtools:q} quickcheck -v {output.bam:q} 2>> {log:q}
        {params.samtools:q} flagstat -@ {threads} {output.bam:q} > {output.flagstat:q}
        """


rule audit_leading_candidate_reads:
    input:
        bam=f"{wgs_root}/alignment/PROBAND01.markdup.bam",
        bai=f"{wgs_root}/alignment/PROBAND01.markdup.bam.bai",
        ranked=f"{private_root}/ranking/mva_prior.ranked.json",
        code="src/mva_hackathon/read_evidence.py",
    output:
        json=f"{wgs_root}/qc/leading_candidate_reads.json",
        review=f"{wgs_root}/qc/leading_candidate_reads.review.md",
    log:
        f"{private_root}/logs/leading_candidate_reads.log",
    params:
        min_mapping_quality=read_audit["min_mapping_quality"],
        min_base_quality=read_audit["min_base_quality"],
        min_phase_support_fragments=read_audit["min_phase_support_fragments"],
    shell:
        """
        mva candidate-read-audit {input.ranked:q} --bam {input.bam:q} \
          --json-output {output.json:q} --review-output {output.review:q} \
          --limit 3 --min-mapping-quality {params.min_mapping_quality} \
          --min-base-quality {params.min_base_quality} \
          --min-phase-support-fragments {params.min_phase_support_fragments} > {log:q}
        """


rule audit_agnostic_candidate_reads:
    input:
        bam=f"{wgs_root}/alignment/PROBAND01.markdup.bam",
        bai=f"{wgs_root}/alignment/PROBAND01.markdup.bam.bai",
        ranked=f"{private_root}/ranking/agnostic_coding.ranked.json",
        code="src/mva_hackathon/read_evidence.py",
    output:
        json=f"{wgs_root}/qc/agnostic_top30_reads.json",
        review=f"{wgs_root}/qc/agnostic_top30_reads.review.md",
    log:
        f"{private_root}/logs/agnostic_top30_reads.log",
    params:
        min_mapping_quality=read_audit["min_mapping_quality"],
        min_base_quality=read_audit["min_base_quality"],
        min_phase_support_fragments=read_audit["min_phase_support_fragments"],
    shell:
        """
        mva candidate-read-audit {input.ranked:q} --bam {input.bam:q} \
          --json-output {output.json:q} --review-output {output.review:q} \
          --limit 30 --min-mapping-quality {params.min_mapping_quality} \
          --min-base-quality {params.min_base_quality} \
          --min-phase-support-fragments {params.min_phase_support_fragments} > {log:q}
        """


rule deepvariant_cpu:
    input:
        reference=wgs_reference,
        fai=f"{wgs_reference}.fai",
        bam=f"{wgs_root}/alignment/PROBAND01.markdup.bam",
        bai=f"{wgs_root}/alignment/PROBAND01.markdup.bam.bai",
    output:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.vcf.gz.tbi",
    threads: config["threads"]["deepvariant"]
    resources:
        mem_mb=80000,
    log:
        f"{private_root}/logs/deepvariant.log",
    params:
        image=deepvariant_image,
        reference_dir=str(Path(wgs_reference).resolve().parent),
        reference_name=Path(wgs_reference).name,
        alignment_dir=str(Path(f"{wgs_root}/alignment").resolve()),
        bam_name="PROBAND01.markdup.bam",
        output_dir=str(Path(f"{wgs_root}/deepvariant").resolve()),
        output_name="PROBAND01.vcf.gz",
        user=f"{os.getuid()}:{os.getgid()}",
    shell:
        """
        mkdir -p {params.output_dir:q}
        docker run --rm --network none --cpus {threads} \
          --user {params.user:q} \
          --volume {params.reference_dir:q}:/reference:ro \
          --volume {params.alignment_dir:q}:/reads:ro \
          --volume {params.output_dir:q}:/output:rw \
          {params.image:q} /opt/deepvariant/bin/run_deepvariant --model_type WGS \
          --ref /reference/{params.reference_name:q} --reads /reads/{params.bam_name:q} \
          --output_vcf /output/{params.output_name:q} \
          --intermediate_results_dir /output/intermediate --num_shards {threads} \
          --vcf_stats_report=true >> {log:q} 2>&1
        """


rule normalize_deepvariant:
    input:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.vcf.gz.tbi",
        reference=wgs_reference,
        fai=f"{wgs_reference}.fai",
    output:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz.tbi",
    threads: config["threads"]["bcftools"]
    log:
        f"{private_root}/logs/deepvariant_normalize.log",
    params:
        bcftools=str(wgs_bin / "bcftools"),
    shell:
        """
        {params.bcftools:q} norm --threads {threads} --multiallelics -any \
          --fasta-ref {input.reference:q} --output-type z --output {output.vcf:q} \
          {input.vcf:q} 2> {log:q}
        {params.bcftools:q} index --force --threads {threads} --tbi \
          {output.vcf:q} 2>> {log:q}
        """


rule split_deepvariant_phase_contig:
    input:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz.tbi",
    output:
        vcf=temp(f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.input.vcf.gz"),
        tbi=temp(f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.input.vcf.gz.tbi"),
    threads: 1
    priority:
        lambda wildcards: 100 if wildcards.phase_contig == "15" else 0
    log:
        f"{private_root}/logs/whatshap_split_{{phase_contig}}.log",
    params:
        bcftools=str(phasing_bin / "bcftools"),
        output_dir=f"{wgs_root}/deepvariant/phasing",
    wildcard_constraints:
        phase_contig=r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)",
    shell:
        """
        mkdir -p {params.output_dir:q}
        {params.bcftools:q} view --regions {wildcards.phase_contig:q} \
          --output-type z --output {output.vcf:q} {input.vcf:q} 2> {log:q}
        {params.bcftools:q} index --force --tbi --output {output.tbi:q} \
          {output.vcf:q} 2>> {log:q}
        """


rule phase_deepvariant_contig:
    input:
        vcf=f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.input.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.input.vcf.gz.tbi",
        bam=f"{wgs_root}/alignment/PROBAND01.markdup.bam",
        bai=f"{wgs_root}/alignment/PROBAND01.markdup.bam.bai",
        reference=wgs_reference,
        fai=f"{wgs_reference}.fai",
    output:
        vcf=temp(f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.phased.vcf.gz"),
        tbi=temp(f"{wgs_root}/deepvariant/phasing/{{phase_contig}}.phased.vcf.gz.tbi"),
    threads: 1
    priority:
        lambda wildcards: 100 if wildcards.phase_contig == "15" else 0
    resources:
        mem_mb=12000,
    log:
        f"{private_root}/logs/whatshap_{{phase_contig}}.log",
    params:
        whatshap=str(phasing_bin / "whatshap"),
        bcftools=str(phasing_bin / "bcftools"),
    wildcard_constraints:
        phase_contig=r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)",
    shell:
        """
        if [ "$({params.bcftools:q} index --nrecords {input.vcf:q})" -eq 0 ]; then
          cp {input.vcf:q} {output.vcf:q}
          : > {log:q}
        else
          {params.whatshap:q} phase --reference {input.reference:q} --sample PROBAND01 \
            --output {output.vcf:q} {input.vcf:q} {input.bam:q} > {log:q} 2>&1
        fi
        {params.bcftools:q} index --force --tbi --output {output.tbi:q} \
          {output.vcf:q} 2>> {log:q}
        """


rule retain_unphased_nonprimary_calls:
    input:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.normalized.vcf.gz.tbi",
    output:
        vcf=temp(f"{wgs_root}/deepvariant/phasing/nonprimary.unphased.vcf.gz"),
        tbi=temp(f"{wgs_root}/deepvariant/phasing/nonprimary.unphased.vcf.gz.tbi"),
    threads: 1
    log:
        f"{private_root}/logs/whatshap_nonprimary.log",
    params:
        bcftools=str(phasing_bin / "bcftools"),
        excluded_primary="^" + ",".join(phase_contigs),
        output_dir=f"{wgs_root}/deepvariant/phasing",
    shell:
        """
        mkdir -p {params.output_dir:q}
        {params.bcftools:q} view --targets {params.excluded_primary:q} \
          --output-type z --output {output.vcf:q} {input.vcf:q} 2> {log:q}
        {params.bcftools:q} index --force --tbi --output {output.tbi:q} \
          {output.vcf:q} 2>> {log:q}
        """


rule phase_deepvariant:
    input:
        phased=phase_shard_vcfs,
        phased_indexes=phase_shard_tbis,
        nonprimary=f"{wgs_root}/deepvariant/phasing/nonprimary.unphased.vcf.gz",
        nonprimary_tbi=f"{wgs_root}/deepvariant/phasing/nonprimary.unphased.vcf.gz.tbi",
    output:
        vcf=f"{wgs_root}/deepvariant/PROBAND01.phased.vcf.gz",
        tbi=f"{wgs_root}/deepvariant/PROBAND01.phased.vcf.gz.tbi",
        stats=f"{wgs_root}/qc/PROBAND01.whatshap.tsv",
    threads: config["threads"]["bcftools"]
    log:
        f"{private_root}/logs/whatshap.log",
    params:
        whatshap=str(phasing_bin / "whatshap"),
        bcftools=str(phasing_bin / "bcftools"),
    shell:
        """
        {params.bcftools:q} concat --threads {threads} --output-type z \
          --output {output.vcf:q} {input.phased:q} {input.nonprimary:q} 2> {log:q}
        {params.bcftools:q} index --force --tbi --output {output.tbi:q} \
          {output.vcf:q} 2>> {log:q}
        {params.whatshap:q} stats --sample PROBAND01 --tsv {output.stats:q} \
          {output.vcf:q} >> {log:q} 2>&1
        """
