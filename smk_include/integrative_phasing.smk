
localrules: master_integrative_phasing,
            write_breakpointr_config_file,
            write_strandphaser_config_file,
            write_phased_vcf_splits_fofn


PATH_INTEGRATIVE_PHASING = '{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}'


rule master_integrative_phasing:
    input:
        []


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
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        strandseq_reads = 'input/fastq/{sseq_reads}.fofn',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/breakpointr.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/breakpointr.input',
    threads: 1
    params:
        bp_cpu = config['num_cpu_high']
    run:
        import os

        outfolder = os.path.dirname(list(input.bam)[0])
        assert os.path.isdir(outfolder), 'write_breakpointr_config_file: alignment output folder does not exist: {}'.format(outfolder)

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
        wc_reg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/{reference}.WCregions.txt',
        cfg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/breakpointR.config',
        rdme = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/README.txt',
        bps = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/breakpoints'),
        bwf = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/browserfiles'),
        data = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/data'),
        plots = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/plots')
    log:
        'log/output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/breakpointr.log'
    benchmark:
        os.path.join('run/output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}',
                     'breakpointr.t{}.rsrc'.format(config['num_cpu_high'])
                     )
    conda:
        '../environment/conda/conda_rscript.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(49152    * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.rdme),
        input_dir = lambda wildcards, input: load_fofn_file(input),
        version = config['git_commit_breakpointr'],
        script_exec = lambda wildcards: find_script_path('run_breakpointr.R')
    shell:
        '{params.script_exec} {params.input_dir} {input.cfg} {params.output_dir} {threads} {output.wc_reg} {params.version} &> {log}'


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
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        strandseq_reads = 'input/fastq/{sseq_reads}.fofn',
        bam = collect_strandseq_alignments  # from module: aux_utilities
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.input',
    threads: 1
    params:
        sp_cpu = config['num_cpu_high']
    run:
        import os

        outfolder = os.path.dirname(list(input.bam)[0])
        assert os.path.isdir(outfolder), 'write_strandphaser_config_file: alignment output folder does not exist: {}'.format(outfolder)

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
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.config',
        fofn = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.input',
        wc_regions = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/{reference}.WCregions.txt',
        variant_calls = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf'
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
        os.path.join('run/output/integrative_phasing/processing/strandphaser',
                     PATH_INTEGRATIVE_PHASING + '.phased.t{}.rsrc'.format(config['num_cpu_high'])
                     )
    conda:
        '../environment/conda/conda_rscript.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(49152 * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        input_dir = lambda wildcards, input: load_fofn_file(input),
        output_dir = lambda wildcards, output: os.path.dirname(output.cfg),
        individual = lambda wildcards: wildcards.sseq_reads.split('_')[0],
        version = config['git_commit_strandphaser'],
        script_exec = lambda wildcards: find_script_path('run_strandphaser.R')
    shell:
        '{params.script_exec} {params.input_dir} {input.cfg} '
            ' {input.variant_calls} {input.wc_regions} '
            ' {params.output_dir} {params.individual} {params.version} &> {log.stp}'


rule write_strandphaser_split_vcf_fofn:
    """
    2021-05-22 - PGAS v13 (needs to be checked if Strand-seq R tools are updated):
    w/o checkpoints, this rule is critical for revealing potential problems for a sample.
    If breakpointR does not identify WC regions for a (sequence) cluster, StrandPhaseR
    won't produce a phased VCF for the respective cluster. This will raise an exception
    in this rule (as intended), and requires the additional check that the "missing" VCF
    is a result of the breakpointR output, and not indicative of an error further upstream in PGAS.

    NB: "reference" here refers to a custom (assembled) reference (clustered squashed assembly)
    """
    input:
        wc_regions = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/{reference}.WCregions.txt',
        vcf_dir = rules.run_strandphaser.output.vcf_dir
    output:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.fofn'
    log:
        'log/output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.write-fofn.log'
    run:
        import os
        import sys

        with open(log[0], 'w') as logfile:
            _ = logfile.write('write_strandphaser_split_vcf_fofn >> processing {}\n'.format(input.vcf_dir))
            # Sanity check: there must be one VCF file per cluster
            sample_name = wildcards.sseq_reads.split('_')[0]
            num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads)
            _ = logfile.write('Estimated number of sequence clusters: {}\n'.format(num_clusters))

            input_dir = input.vcf_dir
            assert os.path.isdir(input_dir), 'write_strandphaser_split_vcf_fofn >> no VCF input folder: {}\n'.format(input_dir)
            input_vcfs = sorted([f for f in os.listdir(input_dir) if f.endswith('.vcf')])
            num_vcf = len(input_vcfs)
            _ = logfile.write('Collected N StrandPhaseR output VCF files: {}\n'.format(num_vcf))
            if num_vcf == 0:
                raise RuntimeError('No StrandPhaseR-phased VCFs in /VCFfiles output dir. '
                                'Cannot create fofn file: {}'.format(output.fofn))

            elif num_vcf != num_clusters:
                _ = logfile.write('Potential error - checking breakpointR output\n')
                brkp_clusters = set()
                with open(input.wc_regions, 'r') as table:
                    for line in table:
                        if not line.strip():
                            continue
                        brkp_clusters.add(line.split()[0])
                _ = logfile.write('breakpointR called W/C-only for N sequence clusters: {}\n'.format(len(brkp_clusters)))
                if len(brkp_clusters) == num_vcf:
                    missing_clusters = set()
                    for i in range(1, num_clusters + 1):
                        cluster_id = 'cluster' + str(i)
                        if cluster_id not in brkp_clusters:
                            missing_clusters.add(cluster_id)
                    missing_clusters = sorted(missing_clusters)
                    if WARN:
                        _ = sys.stderr.write('WARNING: sequence cluster(s) [{}] w/o W/C-only regions for sample: '
                                            '{} / {}\n'.format(', '.join(missing_clusters), wildcards.reference, wildcards.sseq_reads))
                    _ = logfile.write('Sequence cluster(s) w/o W/C-only regions as per breakpointR output: {}\n'.format(', '.join(missing_clusters)))
                else:
                    raise RuntimeError('Mismatch between StrandPhaseR output VCF and number of SaaRclusters: '
                                        '{} VCF files vs {} sequence clusters'.format(num_vcf, num_clusters))
            else:
                pass
            
            with open(output.fofn, 'w') as fofn:
                for file_name in input_vcfs:
                    file_path = os.path.join(input_dir, file_name)
                    _ = fofn.write(file_path + '\n')
            
            _ = logfile.write('Output fofn produced\n')


rule merge_strandphaser_output:
    input:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.fofn'
    output:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.concat-phased.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    shell:
        'bcftools concat -f {input.fofn} --output-type v --output {output} &> {log}'


rule run_integrative_phasing:
    """
    Why this complicated command line?
    One of those "why-is-this-a-problem" situations; for some reason, Snakemake
    does not recognize the output of the rule "run_strandphaser"
    as the producer for the individual VCF splits. So, merge, and split again...
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq',
        spr_phased = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.vcf'
    output:
        vcf = 'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.vcf'
    log:
        'log/output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.log'
    benchmark:
        'run/output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.{sequence}.phased.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'whatshap --debug phase --chromosome {wildcards.sequence} --reference {input.fasta} '
            ' {input.vcf} {input.bam} {input.spr_phased} 2> {log} '
            ' | egrep "^(#|{wildcards.sequence}\s)" > {output}'


def intphase_collect_phased_vcf_split_files(wildcards, glob_collect=True, caller='snakemake'):
    """
    """
    import os

    source_path = os.path.join('output/integrative_phasing/processing/whatshap',
                               PATH_INTEGRATIVE_PHASING,
                               '{hap_reads}.{sequence}.phased.vcf')

    if glob_collect:
        import glob
        pattern = source_path.replace('.{sequence}.', '.*.')
        pattern = pattern.format(**dict(wildcards))
        vcf_files = glob.glob(pattern)

        if not vcf_files:
            if caller == 'snakemake':
                sample_name = wildcards.sseq_reads.split('_')[0]
                num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads)
                tmp = dict(wildcards)
                vcf_files = []
                for i in range(1, num_clusters + 1):
                    tmp['sequence'] = 'cluster' + str(i)
                    vcf_files.append(source_path.format(**tmp))
            else:
                raise RuntimeError('intphase_collect_phased_vcf_split_files: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        folder_path = 'output/reference_assembly/clustered/' + wildcards.sseq_reads
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=folder_path, reference=wildcards.reference).output[0]

        checkpoint_wildcards = glob_wildcards(
            os.path.join(seq_output_dir, '{sequence}.seq')
            )

        vcf_files = expand(
            source_path,
            var_caller=wildcards.var_caller,
            reference=wildcards.reference,
            gq=wildcards.gq,
            qual=wildcards.qual,
            vc_reads=wildcards.vc_reads,
            sseq_reads=wildcards.sseq_reads,
            hap_reads=wildcards.hap_reads,
            sequence=checkpoint_wildcards.sequence
            )

    assert vcf_files, 'intphase_collect_phased_vcf_split_files >> no VCF files returned: {}'.format(wildcards)
    return sorted(vcf_files)


rule write_phased_vcf_splits_fofn:
    """
    2021-05-22 / same as above for StrandPhaseR phased VCF output
    However, here can simplify things assuming that the check above would catch
    any error / prevent pipeline from proceeding.
    """
    input:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_INTEGRATIVE_PHASING + '.spr-phased.fofn',
        vcf_splits = intphase_collect_phased_vcf_split_files
    output:
        fofn = 'output/integrative_phasing/processing/whatshap/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.fofn'
    run:
        import os

        # Sanity check: there must be one VCF file per cluster
        sample_name = wildcards.sseq_reads.split('_')[0]
        num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads)

        num_vcf = len(input.vcf_splits)

        num_spr_vcf = 0
        with open(input.fofn, 'r') as fofn:
            num_spr_vcf = len([line for line in fofn if line.strip()])
        assert num_spr_vcf > 0, 'write_phased_vcf_splits_fofn >> number of StrandPhaseR VCF splits read from fofn is zero: {}'.format(input.fofn)

        if num_vcf == 0:
            raise RuntimeError('write_phased_vcf_splits_fofn >> zero phased VCF split files: {}'.format(wildcards))
        elif num_vcf != num_clusters and num_vcf != num_spr_vcf:
            raise RuntimeError('write_phased_vcf_splits_fofn >> mismatch between expected ({}) and received ({}) phased VCF split files: {}'.format(num_clusters, num_vcf, wildcards))
        else:
            with open(output.fofn, 'w') as dump:
                for file_path in sorted(input.vcf_splits):
                    _ = dump.write(file_path + '\n')


rule strandseq_dga_merge_sequence_phased_vcf_files:
    input:
        fofn = rules.write_phased_vcf_splits_fofn.output.fofn
    output:
        'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.concat.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
             mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
             mem_total_mb = lambda wildcards, attempt: 4096 * attempt
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
    conda:
        '../environment/conda/conda_biotools.yml'
    priority: 200
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
        'bcftools stats {input.vcf} > {output.bcf_stats} 2> {log}'
            ' && '
            'whatshap stats --tsv {output.spr_stats_tsv} {input.vcf} > {output.spr_stats_txt} 2>> {log}'


rule compute_whatshap_phased_vcf_stats:
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz',
        idx = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz.tbi'
    output:
        bcf_stats = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.stats',
        wh_stats_tsv = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.tsv',
        wh_stats_txt = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.txt'
    log:
        'log/output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.stats.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    priority: 200
    resources:
             mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
             mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
        'bcftools stats {input.vcf} > {output.bcf_stats} 2> {log}'
            ' && '
            'whatshap stats --tsv {output.wh_stats_tsv} {input.vcf} > {output.wh_stats_txt} 2>> {log}'


