
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'run_alignments.smk'

localrules: master_variant_calling, \
            merge_sequence_vcf_files

rule master_variant_calling:
    input:


rule compute_position_coverage:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        # require FASTA index here because bedtools fails otherwise
        ref_idx = 'output/reference_assembly/{folder_path}/{reference}.fasta.fai',
        read_ref_aln = 'output/alignments/reads_to_reference/{folder_path}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{folder_path}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai'
    output:
        'output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}.pos-cov.txt.gz'
    benchmark:
        'run/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}.pos-cov.rsrc'
    threads: 2
    shell:
        'bedtools genomecov -d -ibam {input.read_ref_aln} | gzip > {output}'


rule compute_uniform_coverage_regions:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        seq_info = 'output/reference_assembly/{folder_path}/{reference}/sequences/{sequence}.seq',
        pos_cov = rules.compute_position_coverage.output[0]
    output:
        'output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.regions'
    log:
        'log/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.rsrc'
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
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        ref_idx = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta.fai',
        ref_regions = rules.compute_uniform_coverage_regions.output[0]
    output:
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.rsrc'
    params:
        timeout = config['freebayes_timeout_sec'],
        script_dir = config['script_dir']
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 12288,
        mem_total_mb = 393216  # memory consumption varies wildly!!! Replacing FreeBayes at some point would be nice...
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
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        ref_idx = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta.fai',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
    output:
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/00-raw/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/00-raw/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/00-raw/splits/{vc_reads}.{sequence}.rsrc'
    params:
        individual = lambda wildcards: wildcards.vc_reads.split('_')[0]
    resources:
        mem_per_cpu_mb = 16384,
        mem_total_mb = 16384
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
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/00-raw/splits/{vc_reads}.{sequence}.vcf'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.log'
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/norm_longshot_output.py --debug --input-vcf {input} --output-vcf {output} &> {log}'


rule filter_variant_calls_quality_biallelic_snps:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/10-norm/splits/{vc_reads}.{sequence}.vcf.bgz.tbi',
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta'
    output:
        temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf')
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.log'
    shell:
        'bcftools filter --include "QUAL>={wildcards.qual}" {input.vcf} | ' \
            'bcftools view -c 1 --types snps -m 2 -M 2 | ' \
            'bcftools norm --fasta-ref {input.reference} --output {output} --output-type v &> {log}'


rule whatshap_regenotype_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf',
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta'
    output:
       'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.rsrc'
    shell:
        'whatshap genotype --chromosome {wildcards.sequence} --reference {input.reference} --output {output} {input.vcf} {input.read_ref_aln} &> {log}'


rule filter_retyped_by_genotype_quality:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.vcf'
    output:
        temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.vcf')
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.log'
    shell:
        'bcftools filter --include \'GQ>={wildcards.gq}\' --output-type v --output {output} {input.vcf} &> {log}'


rule extract_heterozygous_variants:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf_original = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf',
        vcf_retyped = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.vcf'
    output:
        vcf_original = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf'),
        vcf_retyped = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf'),
    shell:
        'bcftools view --genotype het --output-type v --output-file {output.vcf_original} {input.vcf_original} ' \
        ' && ' \
        'bcftools view --genotype het --output-type v --output-file {output.vcf_retyped} {input.vcf_retyped} '



rule intersect_original_retyped_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        original_vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf.bgz',
        original_tbi = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf.bgz.tbi',
        retyped_vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf.bgz',
        retyped_tbi = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf.bgz.tbi',
    output:
        uniq_original = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0000.vcf'),
        uniq_retyped = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0001.vcf'),
        shared_original = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0002.vcf'),
        shared_retyped = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0003.vcf'),
        desc = temp('output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/README.txt'),
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/{vc_reads}.{sequence}.isect.log'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.desc)
    shell:
        'bcftools isec -p {params.outdir} {input.original_vcf} {input.retyped_vcf} &> {log}'


def collect_final_vcf_splits(wildcards):
    """
    """
    reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                    reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    vcf_files = expand(
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0002.vcf',
        var_caller=wildcards.var_caller,
        reference=wildcards.reference,
        sts_reads=wildcards.sts_reads,
        qual=wildcards.qual,
        gq=wildcards.gq,
        vc_reads=wildcards.vc_reads,
        sequence=checkpoint_wildcards.sequence
        )

    return vcf_files


rule write_final_vcf_splits:
    input:
        vcf_splits = collect_final_vcf_splits
    output:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.fofn'
    run:
        vcf_files = collect_final_vcf_splits(wildcards)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(vcf_files):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for final VCF split merge, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_final_vcf_splits:
    input:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.fofn'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf'
    params:
        vcf_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'bcftools concat --output {output} --output-type v {params.vcf_files}'


rule compute_final_vcf_stats:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf.bgz',
        idx = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf.bgz.tbi',
    output:
        stats = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/stats/{vc_reads}.snps.QUAL{qual}.GQ{gq}.vcf.stats'
    shell:
        'bcftools stats {input.vcf} > {output.stats}'


# Same again for convenience on intermediate set of SNPs

def collect_intermediate_vcf_splits(wildcards):
    """
    """
    reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                    reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    vcf_files = expand(
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/temp/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf',
        var_caller=wildcards.var_caller,
        reference=wildcards.reference,
        sts_reads=wildcards.sts_reads,
        qual=wildcards.qual,
        vc_reads=wildcards.vc_reads,
        sequence=checkpoint_wildcards.sequence
        )

    return vcf_files


rule write_intermediate_vcf_splits:
    input:
        vcf_splits = collect_intermediate_vcf_splits
    output:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}/{vc_reads}.snps.fofn'
    run:
        vcf_files = collect_intermediate_vcf_splits(wildcards)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(vcf_files):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for intermediate VCF split merge, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_intermediate_vcf_splits:
    input:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}/{vc_reads}.snps.fofn'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}/{vc_reads}.snps.vcf'
    params:
        vcf_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'bcftools concat --output {output} --output-type v {params.vcf_files}'


rule compute_intermediate_vcf_stats:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}/{vc_reads}.snps.vcf.bgz',
        idx = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}/{vc_reads}.snps.vcf.bgz.tbi',
    output:
        stats = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/stats/{vc_reads}.snps.QUAL{qual}.vcf.stats'
    shell:
        'bcftools stats {input.vcf} > {output.stats}'


# Deprecated rule...
#rule depth_filter_intersected_variant_calls:
#    """
#    vc_reads = FASTQ file used for variant calling relative to reference
#    """
#    input:
#        'output/variant_calls/{var_caller}/{reference}/split_by_seq/intersect_original_retyped-GQ_{gq}/{vc_reads}.{sequence}.isect.vcf'
#    output:
#        'output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.vcf'
#    log:
#        'log/output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.log'
#    benchmark:
#        'run/output/variant_calls/{var_caller}/{reference}/split_by_seq/final_GQ{gq}_DP{dp}/{vc_reads}.{sequence}.final.rsrc'
#    shell:
#        # TODO: find out why longshot has depth not as FORMAT field
#        "bcftools filter --exclude \'INFO/DP>{wildcards.dp}\' --output-type v --output {output} {input} &> {log}"

