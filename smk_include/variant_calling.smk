
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'run_alignments.smk'

localrules: master_variant_calling

rule master_variant_calling:
    input:


checkpoint create_assembly_sequence_files:
    input:
        'references/assemblies/{reference}.fasta.fai'
    output:
        directory('output/support_files/assembly_sequences/{reference}')
    run:
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                output_path = os.path.join(output_dir, seq_name + '.seq')
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule compute_pos_coverage:
    input:
        ref_idx = 'references/assemblies/{reference}.fasta.fai',
        read_ref_aln = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai'
    output:
        'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}.pos-cov.txt.gz'
    benchmark:
        'run/output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}.pos-cov.rsrc'
    threads: 2
    shell:
        'bedtools genomecov -d -ibam {input.read_ref_aln} | gzip > {output}'


rule compute_uniform_coverage_regions:
    input:
        seq_info = 'output/support_files/assembly_sequences/{reference}/{sequence}.seq',
        pos_cov = 'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}.pos-cov.txt.gz'
    output:
        'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}/{sequence}.unicov.regions'
    log:
        'log/output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}/{sequence}.unicov.log'
    benchmark:
        'run/output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}/{sequence}.unicov.rsrc'
    params:
        num_regions = 128,
        script_dir = config['script_dir']
    shell:
        'zgrep "{wildcards.sequence}\s" {input.pos_cov} | {params.script_dir}/np_cov_to_regions.py --debug ' \
            ' --seq-info {input.seq_info} --num-regions {params.num_regions} --output {output}' \
            ' &> {log}'


def collect_unicov_region_files(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(**wildcards).output[0]
    assert os.path.isdir(seq_output_dir), 'Checkpoint updated, but no output folder at path: {}'.format(seq_output_dir)

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    unicov_region_files = expand(
        'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}/{sequence}.unicov.regions',
        readset=wildcards.readset,
        reference=wildcards.reference,
        sequence=checkpoint_wildcards.sequence
        )
    return unicov_region_files


rule merge_uncov_region_files:
    input:
        regions = collect_unicov_region_files
    output:
        'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}.unicov.regions'
    run:
        # NB: could be too many files to merge directly on the shell
        # due to character limit - the following always works

        with open(output[0], 'w') as dump:
            for regfile in input:
                with open(regfile, 'r') as regions:
                    _ = dump.write(regions.read())
    # rule end


rule call_variants_freebayes_parallel:
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta',
        ref_idx = 'references/assemblies/{reference}.fasta.fai',
        ref_regions = 'output/support_files/freebayes_unicov_regions/{readset}_map-to_{reference}.unicov.regions'
    output:
        'output/variant_calls/freebayes/{reference}/norm/{readset}.vcf'
    log:
        'log/output/variant_calls/freebayes/{reference}/norm/{readset}.parallel.log'
    benchmark:
        'run/output/variant_calls/freebayes/{reference}/norm/{readset}.rsrc'
    wildcard_constraints:
        reference = '[\w\-]+'
    params:
        timeout = 3600,  # timeout in seconds
        script_dir = config['script_dir']
    threads: 128
    shell:
        '{params.script_dir}/fb-parallel-timeout.sh {input.ref_regions} {threads} {params.timeout} {log}' \
            ' --use-best-n-alleles 4 -f {input.reference} {input.read_ref_aln} > {output}'


rule call_variants_longshot:
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta',
        ref_idx = 'references/assemblies/{reference}.fasta.fai'
    output:
        'output/variant_calls/longshot/{reference}/raw/{readset}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/raw/{readset}.log'
    benchmark:
        'run/output/variant_calls/longshot/{reference}/raw/{readset}.rsrc'
    params:
        individual = lambda wildcards: wildcards.readset.split('_')[0]
    shell:
        'longshot --no_haps --bam {input.read_ref_aln} --ref {input.reference} --sample_id {params.individual} --out {output} &> {log}'


rule normalize_longshot_vcf:
    input:
        'output/variant_calls/longshot/{reference}/raw/{readset}.vcf'
    output:
        'output/variant_calls/longshot/{reference}/norm/{readset}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/norm/{readset}.log'
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/norm_longshot_output.py --debug --input-vcf {input} --output-vcf {output} &> {log}'


rule quality_filter_variant_calls:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/norm/{readset}.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/norm/{readset}.vcf.bgz.tbi',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.rsrc'
    run:
        exec = "bcftools filter"
        exec += " --include \'QUAL>=10\' {input.vcf}"
        exec += " | bcftools view -c 1 --types snps -m 2 -M 2"
        exec += " | bcftools norm --fasta-ref {input.reference} --output {output} --output-type v"
        exec += " &> {log}"
        shell(exec)


rule whatshap_regenotype_variant_calls:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.vcf',
        read_ref_aln = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/variant_calls/{var_caller}/{reference}/regenotype_whatshap/{readset}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/regenotype_whatshap/{readset}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/regenotype_whatshap/{readset}.rsrc'
    shell:
        'whatshap genotype --reference {input.reference} --output {output} {input.vcf} {input.read_ref_aln} &> {log}'


rule genotype_quality_filter_retyped:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/regenotype_whatshap/{readset}.vcf',
    output:
        'output/variant_calls/{var_caller}/{reference}/filter_regenotyped_GQ_{gq}/{readset}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/filter_regenotyped_GQ_{gq}/{readset}.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/filter_regenotyped_GQ_{gq}/{readset}.rsrc'
    shell:
        'bcftools filter --include \'GQ>={wildcards.gq}\' --output-type v --output {output} {input.vcf} &> {log}'


rule intersect_original_retyped_variant_calls:
    input:
        original_vcf = 'output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.het-only.vcf.bgz',
        original_tbi = 'output/variant_calls/{var_caller}/{reference}/filter_qual_type/{readset}.het-only.vcf.bgz.tbi',
        retyped_vcf = 'output/variant_calls/{var_caller}/{reference}/filter_regenotyped_GQ_{gq}/{readset}.het-only.vcf.bgz',
        retyped_tbi = 'output/variant_calls/{var_caller}/{reference}/filter_regenotyped_GQ_{gq}/{readset}.het-only.vcf.bgz.tbi'
    output:
        uniq_original = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}/0000.vcf',
        uniq_retyped = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}/0001.vcf',
        shared_original = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}/0002.vcf',
        shared_retyped = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}/0003.vcf',
        desc = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}/README.txt',
        readset_copy = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}.isect.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}.isect.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}.isect.rsrc'
    params:
        outdir = 'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}'
    shell:
        'bcftools isec -p {params.outdir} {input.original_vcf} {input.retyped_vcf} &> {log} && cp {output.shared_original} {output.readset_copy}'


rule depth_filter_intersected_variant_calls:
    input:
        'output/variant_calls/{var_caller}/{reference}/intersect_original_retyped-GQ_{gq}/{readset}.isect.vcf'
    output:
        'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{readset}.final.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{readset}.final.log'
    benchmark:
        'run/output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{readset}.final.rsrc'
    shell:
        # TODO: find out why longshot has depth not as FORMAT field
        "bcftools filter --exclude \'INFO/DP>{wildcards.dp}\' --output-type v --output {output} {input} &> {log}"


rule extract_heterozygous_variants:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/filter_{filter_type}/{readset}.vcf',
    output:
        'output/variant_calls/{var_caller}/{reference}/filter_{filter_type}/{readset}.het-only.vcf',
    shell:
        'bcftools view --genotype het --output-type v --output-file {output} {input.vcf}'
