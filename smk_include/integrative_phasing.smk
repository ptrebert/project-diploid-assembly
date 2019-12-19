
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
    conda:
        config['conda_env_rdga']
    resources:
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
    conda:
        config['conda_env_rdga']
    resources:
        mem_total_mb = 3072,
        mem_per_cpu_mb = 3072
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_strandphaser']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_strandphaser.R {params.version} &> {log}'


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
        strandseq_reads = 'input/fastq/strand-seq/{sts_reads}.fofn',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/breakpointr.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/breakpointr.input',
    threads: 1
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
        wc_reg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        cfg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/breakpointR.config',
        rdme = 'output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/README.txt',
        bps = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/breakpoints'),
        bwf = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/browserfiles'),
        data = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/data'),
        plots = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/run/plots')
    log:
        'log/output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/breakpointr.log'
    benchmark:
        'run/output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/breakpointr.rsrc'
    conda:
        config['conda_env_rdga']
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.rdme),
        input_dir = lambda wildcards, input: load_fofn_file(input),
        script_dir = config['script_dir']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(3072 / config['num_cpu_high']),
        mem_total_mb = 3072,
        runtime_hrs = 6
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
        strandseq_reads = 'input/fastq/strand-seq/{sts_reads}.fofn',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/strandphaser.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/strandphaser.input',
    threads: 1
    params:
        sp_cpu = config['num_cpu_high']
    run:
        import os
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
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/strandphaser.config',
        fofn = 'output/integrative_phasing/processing/config_files/{reference}/{sts_reads}/strandphaser.input',
        wc_regions = 'output/integrative_phasing/processing/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        variant_calls = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf'
    output:
        browser = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/browserFiles'),
        data = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/data'),
        phased = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/Phased'),
        maps = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/SingleCellHaps'),
        vcf_dir = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/VCFfiles'),
        cfg = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/StrandPhaseR.config',
    log:
        stp = 'log/output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.phased.log',
    benchmark:
        'run/output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.phased.rsrc'
    conda:
        config['conda_env_rdga']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(8192 / config['num_cpu_high']),
        mem_total_mb = 8192,
        runtime_hrs = 6
    params:
        input_dir = lambda wildcards, input: load_fofn_file(input),
        output_dir = lambda wildcards, output: os.path.dirname(output.cfg),
        individual = lambda wildcards: wildcards.sts_reads.split('_')[0],
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/run_strandphaser.R {params.input_dir} {input.cfg} '
            ' {input.variant_calls} {input.wc_regions} '
            ' {params.output_dir} {params.individual} &> {log.stp}'


rule normalize_strandphaser_output:
    """
    StrandPhaseR does not report unobserved alleles, which creates
    incompatible VCF records for WhatsHap (records with missing genotype
    will be ignored). Since, in our case, we assume all variants to be
    het SNVs, simply complement if missing
    """
    input:
        rules.run_strandphaser.output.vcf_dir
    output:
        norm_dir = directory('output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '/strandphaser_norm'),
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.fofn'
    resources:
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048,
    message:
        "Running potentially deprecated task - after commit #e608407, "
           "StrandPhaseR does no longer omit unobserved alleles in output VCF."
    run:
        import io
        import os
        os.makedirs(output.norm_dir, exist_ok=True)

        input_dir = input[0]
        input_vcfs = [f for f in os.listdir(input_dir) if f.endswith('.vcf')]
        norm_paths = []

        gt_map = {'1|.': '1|0',
                  '0|.': '0|1',
                  '.|1': '0|1',
                  '.|0': '1|0'}

        for vcf_file in input_vcfs:
            input_path = os.path.join(input_dir, vcf_file)
            vcf_buffer = io.StringIO()
            with open(input_path, 'r') as vcf:
                for line in vcf:
                    if not line.startswith('#'):
                        columns = line.strip().split('\t')
                        format_colum = columns[-2]
                        value_column = columns[-1]
                        assert format_colum.startswith('GT:'), \
                            'Unexpected StrandPhaseR output format: {} / {}'.format(input_path, line.strip())
                        values = value_column.split(':')
                        values[0] = gt_map.get(values[0], values[0])
                        value_column = ':'.join(values)
                        columns[-1] = value_column
                        line = '\t'.join(columns) + '\n'
                    vcf_buffer.write(line)
                output_path = os.path.join(output.norm_dir, vcf_file)
                with open(output_path, 'w') as dump:
                    _ = dump.write(vcf_buffer.getvalue())
                norm_paths.append(output_path)

        with open(output.fofn, 'w') as dump:
            _ = dump.write('\n'.join(sorted(norm_paths)))


rule merge_normalized_strandphaser_output:
    input:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.fofn'
    output:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.concat-phased.log'
    shell:
        'bcftools concat -f {input.fofn} --output-type v --output {output} &> {log}'


rule run_integrative_phasing:
    """
    Why this complicated command line?
    One of those "why-is-this-a-problem" situations; for some reason, Snakemake
    does not recognize the output of the rule "normalize_strandphaser_output"
    as the producer for the individual VCF splits. So, merge, and split again...
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/{sts_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
        spr_phased = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf'
    output:
        vcf = 'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.vcf'
    log:
        'log/output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.log'
    benchmark:
        'run/output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.rsrc'
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = 2048
    shell:
        'whatshap --debug phase --chromosome {wildcards.sequence} --reference {input.fasta} '
            ' {input.vcf} {input.bam} {input.spr_phased} 2> {log} '
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
        'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.vcf',
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
        fofn = 'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.fofn'
    run:
        # follow same example as merge strand-seq BAMs in module prepare_custom_references
        vcf_files = intphase_collect_phased_vcf_split_files(wildcards)
        if len(vcf_files) == 0:
            raise RuntimeError('No phased VCF files to merge. Previous job(s) likely failed for: '
                               '{}'.format(output.fofn))

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
        'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.concat.log'
    shell:
        'bcftools concat -f {input.fofn} --output {output} --output-type v &> {log}'


rule compute_strandphaser_phased_vcf_stats:
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf.bgz',
        idx = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf.bgz.tbi'
    output:
        bcf_stats = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.vcf.stats',
        spr_stats_tsv = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.stats.tsv',
        spr_stats_txt = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.stats.txt'
    log:
        'log/output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.stats.log'
    shell:
        'bcftools stats {input.vcf} > {output.bcf_stats} 2> {log}'
            ' && '
            'whatshap stats --tsv {output.spr_stats_tsv} {input.vcf} > {output.spr_stats_txt} 2>> {log}'


rule compute_whatshap_phased_vcf_stats:
    input:
        spr_bcf_stats = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.vcf.stats',
        spr_stats_tsv = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.stats.tsv',
        spr_stats_txt = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.spr-phased.stats.txt',
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz',
        idx = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz.tbi'
    output:
        bcf_stats = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.stats',
        wh_stats_tsv = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.tsv',
        wh_stats_txt = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.txt'
    log:
        'log/output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.log'
    shell:
        'bcftools stats {input.vcf} > {output.bcf_stats} 2> {log}'
            ' && '
            'whatshap stats --tsv {output.wh_stats_tsv} {input.vcf} > {output.wh_stats_txt} 2>> {log}'


