
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_integrative_phasing


PATH_INTEGRATIVE_PHASING = '{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}'


rule master_integrative_phasing:
    input:


rule install_rlib_breakpointr:
    input:
        'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust'])
    output:
         check = touch('output/check_files/R_setup/breakpointr.ok')
    log:
        'log/output/check_files/R_setup/breakpointr.log'
    resources:
        runtime_hrs = 0,
        runtime_min = 30,
        mem_total_mb = 3072,
        mem_per_cpu_mb = 3072
    params:
        script_dir = config['script_dir']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_breakpointr.R &> {log}'


rule install_rlib_strandphaser:
    input:
        rules.install_rlib_breakpointr.output.check
    output:
         check = touch('output/check_files/R_setup/strandphaser_ver-{}.ok'.format(config['git_commit_strandphaser']))
    log:
        'log/output/check_files/R_setup/strandphaser_ver-{}.log'.format(config['git_commit_strandphaser'])
    resources:
        runtime_hrs = 0,
        runtime_min = 30,
        mem_total_mb = 3072,
        mem_per_cpu_mb = 3072
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_strandphaser']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_strandphaser.R {params.version} &> {output}'


rule write_breakpointr_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for breakpointR
    """
    input:
        setup_ok = rules.install_rlib_breakpointr.output.check,
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        strandseq_reads = 'input/fastq/complete/{sts_reads}.fastq.gz',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/breakpointr.config',
        input_dir = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/breakpointr.input',
    threads: 1
    resources:
        runtime_hrs = 0,
        runtime_min = 5
    params:
        bp_cpu = config['num_cpu_high']
    run:
        # following same example as merge strand-seq BAMs in module prepare_custom_references
        bam_files = collect_strandseq_alignments(wildcards)
        outfolder = os.path.dirname(bam_files[0])

        config_rows = [
            '[General]',
            'numCPU = ' + str(params.bp_cpu),  # due to a bug in breakpointr, this value has to be repeated on the CLI
            'reuse.existing.files = FALSE',
            '',
            '[breakpointR]',
            'windowsize = 500000',
            'binMethod = "size"',
            'pairedEndReads = TRUE',
            'pair2frgm = FALSE',
            'min.mapq = 10',
            'filtAlt = TRUE',
            'background = 0.1',
            'minReads = 50'
        ]

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_breakpointr:
    input:
        cfg = rules.write_breakpointr_config_file.output.cfg,
        fofn = rules.write_breakpointr_config_file.output.input_dir
    output:
        wc_reg = 'output/integrative_phasing/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        cfg = 'output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/breakpointR.config',
        rdme = 'output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/README.txt',
        bps = directory('output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/breakpoints'),
        bwf = directory('output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/browserfiles'),
        data = directory('output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/data'),
        plots = directory('output/integrative_phasing/breakpointr/{reference}/{sts_reads}/run/plots')
    log:
        'log/output/integrative_phasing/breakpointr/{reference}/{sts_reads}/breakpointr.log'
    benchmark:
        'run/output/integrative_phasing/breakpointr/{reference}/{sts_reads}/breakpointr.rsrc'
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.rdme),
        input_dir = lambda wildcards, input: load_fofn_file(input),
        script_dir = config['script_dir']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(3072 / config['num_cpu_high']),
        mem_total_mb = 3072,
        runtime_hrs = 5
    shell:
        '{params.script_dir}/run_breakpointr.R {params.input_dir} {input.cfg} {params.output_dir} {threads} {output.wc_reg} &> {log}'


rule write_strandphaser_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for StrandPhaseR
    """
    input:
        setup_ok = rules.install_rlib_strandphaser.output.check,
        reference = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        strandseq_reads = 'input/fastq/complete/{sts_reads}.fastq.gz',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.config',
        input_dir = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.input',
    threads: 1
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    params:
        sp_cpu = config['num_cpu_high']
    run:
        # following same example as merge strand-seq BAMs in module prepare_custom_references
        bam_files = collect_strandseq_alignments(wildcards)
        outfolder = os.path.dirname(bam_files[0])

        config_rows = [
            '[General]',
            'numCPU = ' + str(params.sp_cpu),
            'pairedEndReads = TRUE',
            'min.mapq = 10',
            '',
            '[StrandPhaseR]',
            'min.baseq = 20',
            'num.iterations = 2',
            'translateBases = TRUE',
            'splitPhasedReads = TRUE'
        ]

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_strandphaser:
    input:
        cfg = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.config',
        fofn = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.input',
        wc_regions = 'output/integrative_phasing/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        variant_calls = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf'
    output:
        browser = directory('output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/browserFiles'),
        data = directory('output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/data'),
        phased = directory('output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/Phased'),
        maps = directory('output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/SingleCellHaps'),
        vcf_dir = directory('output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/VCFfiles'),
        cfg = 'output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/run/StrandPhaseR.config',
        vcf = 'output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.phased.vcf'
    log:
        stp = 'log/output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.phased.log',
        bcf = 'log/output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.concat.log',
    benchmark:
        'run/output/integrative_phasing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.phased.rsrc'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = int(8192 / config['num_cpu_high']),
        mem_total_mb = 8192
    params:
        input_dir = lambda wildcards, input: load_fofn_file(input),
        output_dir = lambda wildcards, output: os.path.dirname(output.cfg),
        individual = lambda wildcards: wildcards.sts_reads.split('_')[0],
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/run_strandphaser.R {params.input_dir} {input.cfg} ' \
            ' {input.variant_calls} {input.wc_regions} ' \
            ' {params.output_dir} {params.individual} &> {log.stp}' \
            ' && ' \
            ' bcftools concat --output {output.vcf} --output-type v `ls {output.vcf_dir}/*phased.vcf` ' \
            ' &> {log.bcf}'


rule run_integrative_phasing:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snps.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
        sts_phased = rules.run_strandphaser.output.vcf
    output:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/splits/{hap_reads}.{sequence}.phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/splits/{hap_reads}.{sequence}.phased.log'
    benchmark:
        'run/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/splits/{hap_reads}.{sequence}.phased.rsrc'
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = 2048
    shell:
        'whatshap --debug phase --chromosome {wildcards.sequence} --reference {input.fasta} ' \
            ' {input.vcf} {input.bam} {input.sts_phased} 2> {log} ' \
            ' | egrep "^(#|{wildcards.sequence}\s)" > {output}'


def intphase_collect_phased_vcf_split_files(wildcards):
    """
    """
    folder_path = 'output/reference_assembly/clustered/' + wildcards.sts_reads
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=folder_path, reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    vcf_files = expand(
        'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/splits/{hap_reads}.{sequence}.phased.vcf',
        var_caller=wildcards.var_caller,
        reference=wildcards.reference,
        gq=wildcards.gq,
        qual=wildcards.qual,
        vc_reads=wildcards.vc_reads,
        sts_reads=wildcards.sts_reads,
        hap_reads=wildcards.hap_reads,
        sequence=checkpoint_wildcards.sequence
        )
    return sorted(vcf_files)


rule write_phased_vcf_splits_fofn:
    input:
        splits = intphase_collect_phased_vcf_split_files
    output:
        fofn = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.phased.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 5
    run:
        # follow same example as merge strand-seq BAMs in module prepare_custom_references
        vcf_files = intphase_collect_phased_vcf_split_files(wildcards)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(vcf_files):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for INTPHASE VCF split merge, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule strandseq_dga_merge_sequence_phased_vcf_files:
    input:
        fofn = rules.write_phased_vcf_splits_fofn.output.fofn
    output:
        'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.phased.vcf'
    resources:
        runtime_hrs = 0,
        runtime_min = 30
    params:
        merge_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'bcftools concat --output {output} --output-type v {params.merge_files}'


rule compute_phased_vcf_stats:
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.phased.vcf.bgz',
        idx = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.phased.vcf.bgz.tbi'
    output:
        stats = 'output/statistics/variant_calls/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.snps.phased.vcf.stats'
    resources:
        runtime_hrs = 0,
        runtime_min = 30
    shell:
        'bcftools stats {input.vcf} > {output.stats}'