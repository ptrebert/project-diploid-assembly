
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_integrative_phasing, \
            strandseq_dga_merge_sequence_phased_vcf_files, \
            install_rlib_breakpointr, \
            install_rlib_strandphaser, \
            write_breakpointr_config_file, \
            write_strandphaser_config_file


PATH_INTEGRATIVE_PHASING = '{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}'


rule master_integrative_phasing:
    input:


rule install_rlib_breakpointr:
    input:
        'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust'])
    output:
         check = 'output/check_files/R_setup/breakpointr.ok'
    params:
        script_dir = config['script_dir']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_breakpointr.R &> {output}'


rule install_rlib_strandphaser:
    input:
        rules.install_rlib_breakpointr.output.check
    output:
         check = 'output/check_files/R_setup/strandphaser_ver-{}.ok'.format(config['git_commit_strandphaser'])
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
    params:
        bp_cpu = config['num_cpu_high']
    run:
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

        outfolder = os.path.dirname(input.bam[0])
        assert os.path.isdir(outfolder), 'BreakpointR: invalid output folder / strand-seq alignments: {}'.format(outfolder)
        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_breakpointr:
    input:
        cfg = rules.write_breakpointr_config_file.output.cfg,
        dir = rules.write_breakpointr_config_file.output.input_dir
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
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        input_dir = lambda wildcards, input: open(input.dir).readline().strip(),
        script_dir = config['script_dir']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 128,
        mem_total_mb = config['num_cpu_high'] * 128
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
    params:
        sp_cpu = config['num_cpu_high']
    run:
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

        outfolder = os.path.dirname(input.bam[0])
        assert os.path.isdir(outfolder), 'StrandPhaseR: invalid output folder / strand-seq alignments: {}'.format(outfolder)
        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_strandphaser:
    input:
        cfg = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.config',
        dir = 'output/integrative_phasing/config_files/{reference}/{sts_reads}/strandphaser.input',
        wc_regions = 'output/integrative_phasing/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        variant_calls = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf'
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
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 256,
        mem_total_mb = config['num_cpu_high'] * 256
    params:
        input_dir = lambda wildcards, input: open(input.dir).readline().strip(),
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
        vcf = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz.tbi',
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
        dp=wildcards.dp,
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
    run:
        potential_log = os.path.join('log', output.fofn.replace('.fofn', '.log'))
        os.makedirs(os.path.dirname(potential_log), exist_ok=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(input.splits):
                if not os.path.isfile(file_path):
                    with open(potential_log, 'w') as error_log:
                        _ = error_log.write('Invalid path to merge VCF split file: {}\n'.format(file_path))
                        _ = error_log.write('Input BAMS: {}\n'.format(input.bams))
                        _ = error_log.write('Type: {}\n'.format(type(input.bams)))
                        raise AssertionError('Invalid path to merge VCF split: {}\n'.format(potential_log))
                _ = dump.write(file_path + '\n')


rule strandseq_dga_merge_sequence_phased_vcf_files:
    input:
        fofn = rules.write_phased_vcf_splits_fofn.output.fofn
    output:
        'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.phased.vcf'
    params:
        merge_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'bcftools concat --output {output} --output-type v {params.merge_files}'
