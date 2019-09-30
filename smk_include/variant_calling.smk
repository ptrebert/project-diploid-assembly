
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'run_alignments.smk'

localrules: master_variant_calling

rule master_variant_calling:
    input:


rule compute_position_coverage:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        ref_idx = 'references/assemblies/{reference}.fasta.fai',
        read_ref_aln = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam.bai'
    output:
        'output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}.pos-cov.txt.gz'
    benchmark:
        'run/output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}.pos-cov.rsrc'
    threads: 2
    shell:
        'bedtools genomecov -d -ibam {input.read_ref_aln} | gzip > {output}'


rule compute_uniform_coverage_regions:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        seq_info = 'references/assemblies/{reference}/sequences/{sequence}.seq',
        pos_cov = 'output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}.pos-cov.txt.gz'
    output:
        'output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.regions'
    log:
        'log/output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.log'
    benchmark:
        'run/output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.rsrc'
    wildcard_constraints:
        vc_reads = '[\w\-]+',
        reference = '[\w\-]+',
        sequence = '\w+'
    params:
        num_regions = 128,
        script_dir = config['script_dir']
    shell:
        'zgrep "{wildcards.sequence}\s" {input.pos_cov} | {params.script_dir}/np_cov_to_regions.py --debug ' \
            ' --seq-info {input.seq_info} --num-regions {params.num_regions} --output {output}' \
            ' &> {log}'


rule call_variants_freebayes_parallel:
    """
    NB: This works on separate region files, and in parallel on each file.
        In other words, for a substantial speed up, this needs to be
        distributed across a cluster or, at least, this job should run
        on a large compute node

    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta',
        ref_idx = 'references/assemblies/{reference}.fasta.fai',
        ref_regions = 'output/alignments/reads_to_reference/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.regions'
    output:
        'output/variant_calls/freebayes/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/freebayes/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.parallel.log'
    benchmark:
        'run/output/variant_calls/freebayes/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.rsrc'
    wildcard_constraints:
        vc_reads = '[\w\-]+',
        reference = '[\w\-]+',
        sequence = '\w+'
    params:
        timeout = 2700,  # timeout in seconds
        script_dir = config['script_dir']
    threads: 16
    shell:
        '{params.script_dir}/fb-parallel-timeout.sh {input.ref_regions} {threads} {params.timeout} {log}' \
            ' --use-best-n-alleles 4 -f {input.reference} {input.read_ref_aln} > {output}'


rule call_variants_longshot:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta',
        ref_idx = 'references/assemblies/{reference}.fasta.fai',
        ref_split = 'references/assemblies/{reference}/sequences/{sequence}.seq',
    output:
        'output/variant_calls/longshot/{reference}/split_by_seq/raw/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/split_by_seq/raw/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/longshot/{reference}/split_by_seq/raw/{vc_reads}.{sequence}.rsrc'
    params:
        individual = lambda wildcards: wildcards.vc_reads.split('_')[0]
    shell:
        'longshot --no_haps --bam {input.read_ref_aln} ' \
            ' --ref {input.reference} --region {wildcards.sequence}' \
             ' --sample_id {params.individual} --out {output} &> {log}'


rule normalize_longshot_vcf:
    """
    NB: Currently, longshot does not adhere to the VCF specification
        and uses a "." as decimal mark for the GQ format field. WhatsHap
        cannot handle this. This step simply rounds the values to full
        integers.

    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        'output/variant_calls/longshot/{reference}/split_by_seq/raw/{vc_reads}.{sequence}.vcf'
    output:
        'output/variant_calls/longshot/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.log'
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/norm_longshot_output.py --debug --input-vcf {input} --output-vcf {output} &> {log}'


rule quality_filter_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/norm/{vc_reads}.{sequence}.vcf.bgz.tbi',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.rsrc'
    wildcard_constraints:
        vc_reads = '[\w\-]+',
        reference = '[\w\-]+',
        sequence = '\w+'
    run:
        exec = "bcftools filter"
        exec += " --include \'QUAL>=10\' {input.vcf}"
        exec += " | bcftools view -c 1 --types snps -m 2 -M 2"
        exec += " | bcftools norm --fasta-ref {input.reference} --output {output} --output-type v"
        exec += " &> {log}"
        shell(exec)


rule whatshap_regenotype_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.vcf',
        read_ref_aln = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/regenotype_whatshap/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/regenotype_whatshap/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/regenotype_whatshap/{vc_reads}.{sequence}.rsrc'
    shell:
        'whatshap genotype --chromosome {wildcards.sequence} --reference {input.reference} --output {output} {input.vcf} {input.read_ref_aln} &> {log}'


rule genotype_quality_filter_retyped:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/regenotype_whatshap/{vc_reads}.{sequence}.vcf'
    output:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_regenotyped_GQ_{gq}/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_regenotyped_GQ_{gq}/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_regenotyped_GQ_{gq}/{vc_reads}.{sequence}.rsrc'
    shell:
        'bcftools filter --include \'GQ>={wildcards.gq}\' --output-type v --output {output} {input.vcf} &> {log}'


rule intersect_original_retyped_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        original_vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.het-only.vcf.bgz',
        original_tbi = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_qual_type/{vc_reads}.{sequence}.het-only.vcf.bgz.tbi',
        retyped_vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_regenotyped_GQ_{gq}/{vc_reads}.{sequence}.het-only.vcf.bgz',
        retyped_tbi = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_regenotyped_GQ_{gq}/{vc_reads}.{sequence}.het-only.vcf.bgz.tbi'
    output:
        uniq_original = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}/0000.vcf',
        uniq_retyped = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}/0001.vcf',
        shared_original = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}/0002.vcf',
        shared_retyped = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}/0003.vcf',
        desc = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}/README.txt',
        readset_copy = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}.isect.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}.isect.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}.isect.rsrc'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.desc)
    shell:
        'bcftools isec -p {params.outdir} {input.original_vcf} {input.retyped_vcf} &> {log} && cp {output.shared_original} {output.readset_copy}'


rule depth_filter_intersected_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}.isect.vcf'
    output:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.rsrc'
    shell:
        # TODO: find out why longshot has depth not as FORMAT field
        "bcftools filter --exclude \'INFO/DP>{wildcards.dp}\' --output-type v --output {output} {input} &> {log}"


rule extract_heterozygous_variants:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_{filter_type}/{vc_reads}.{sequence}.vcf',
    output:
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/filter_{filter_type}/{vc_reads}.{sequence}.het-only.vcf',
    wildcard_constraints:
        vc_reads = '[\w\-]+'
    shell:
        'bcftools view --genotype het --output-type v --output-file {output} {input.vcf}'


def collect_sequence_vcf_files(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(**wildcards).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    vcf_files = expand(
        'output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.vcf',
        var_caller=wildcards.var_caller,
        reference=wildcards.reference,
        gq=wildcards.gq,
        dp=wildcards.dp,
        vc_reads=wildcards.vc_reads,
        sequence=checkpoint_wildcards.sequence
        )
    return vcf_files


rule merge_sequence_vcf_files:
    input:
        vcf_files = collect_sequence_vcf_files
    output:
        'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf',
    wildcard_constraints:
        var_caller = '(freebayes|longshot)',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        dp = '[0-9]+',
        gq = '[0-9]+'
    shell:
        'bcftools concat --output {output} --output-type v {input.vcf_files}'