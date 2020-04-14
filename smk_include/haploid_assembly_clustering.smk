
import os

localrules: master_haploid_assembly_clustering


rule master_haploid_assembly_clustering:
    input:
        []


def hac_collect_strandseq_merge_files(wildcards, glob_collect=False):
    """
    """
    source_path = os.path.join(
        'output/alignments/strandseq_to_phased_assembly',
        'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
        '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
        'temp/aln/{individual}_{project}_{platform}-{spec}_{lib_id}_{run_id}.filt.sam.bam')

    individual = wildcards.individual
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    assert individual in wildcards.reference and individual in wildcards.sseq_reads, \
        'Wrong reference / sseq_reads match: {} / {}'.format(wildcards.reference, wildcards.sseq_reads)

    if glob_collect:
        import glob
        pattern = source_path.replace('{run_id}', '*')
        pattern = pattern.replace('{spec}', '*')
        pattern = pattern.format(**dict(wildcards))
        bam_files = glob.glob(pattern)

        if not bam_files:
            raise RuntimeError('hac_collect_strandseq_merge_files: no files collected with pattern {}'.format(pattern))

    else:
        requests_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sseq_reads).output[0]

        search_pattern = '_'.join([individual, project, platform + '-{spec}', lib_id, '{run_id}', '1'])
        search_path = os.path.join(requests_dir, search_pattern + '.request')

        checkpoint_wildcards = glob_wildcards(search_path)

        bam_files = expand(
            source_path,
            zip,
            hap_assm_mode=[wildcards.hap_assm_mode, wildcards.hap_assm_mode],
            var_caller=[wildcards.var_caller, wildcards.var_caller],
            qual=[wildcards.qual, wildcards.qual],
            gq=[wildcards.gq, wildcards.gq],
            reference=[wildcards.reference, wildcards.reference],
            vc_reads=[wildcards.vc_reads, wildcards.vc_reads],
            sseq_reads=[wildcards.sseq_reads, wildcards.sseq_reads],
            pol_reads=[wildcards.pol_reads, wildcards.pol_reads],
            hap_reads=[wildcards.hap_reads, wildcards.hap_reads],
            hap_assembler=[wildcards.hap_assembler, wildcards.hap_assembler],
            hap=[wildcards.hap, wildcards.hap],
            pol_pass=[wildcards.pol_pass, wildcards.pol_pass],
            individual=[individual, individual],
            project=[project, project],
            platform=[platform, platform],
            spec=checkpoint_wildcards.spec,
            lib_id=[lib_id, lib_id],
            run_id=checkpoint_wildcards.run_id)

    return bam_files


rule hac_write_strandseq_merge_fofn:
    """
    """
    input:
        bams = hac_collect_strandseq_merge_files
    output:
        fofn = os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            'temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.fofn'
        )
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES,
        #lib_id = 'P[A-Z0-9]+'
    run:
        import os
        pattern = '[empty]'
        try:
            validate_checkpoint_output(input.bams)
            bam_files = input.bams
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            bam_files = hac_collect_strandseq_merge_files(wildcards, glob_collect=True)

        if len(bam_files) != 2:
            raise RuntimeError('(hac) Missing merge partner for strand-seq BAM files {} / {}: '
                               '{}'.format(output.fofn, pattern, bam_files))

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(bam_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\n(hac) WARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule hac_merge_mono_dinucleotide_fraction:
    """
    This rule is likely quite specific
    for the strand-seq data used in this
    pipeline - point of 'failure' for
    different input data
    """
    input:
        fofn = rules.hac_write_strandseq_merge_fofn.output.fofn
    output:
        temp(
            os.path.join(
                'output/alignments/strandseq_to_phased_assembly',
                'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
                '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
                'temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam'
            )
        )
    log:
        os.path.join(
            'log', 'output/alignments/strandseq_to_phased_assembly',
                'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
                '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
                'temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.log'
        )
    benchmark:
        os.path.join(
            'run', 'output/alignments/strandseq_to_phased_assembly',
                'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
                '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}', 'temp/mrg',
                '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg' + '.t{}.rsrc'.format(config['num_cpu_low'])
        )
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    params:
        merge_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'samtools merge -@ {threads} -O BAM {output} {params.merge_files} &> {log}'


rule hac_link_strandseq_monofraction_samples:
    input:
        os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            'temp/aln/{individual}_{project}_{platform}-{spec}_{lib_id}.filt.sam.bam')
    output:
        os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            'temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam')
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES
    params:
        shell_cmd = lambda wildcards: 'cp' if bool(config.get('force_local_copy', False)) else 'ln -s'
    shell:
        '{params.shell_cmd} {input} {output}'


rule hac_samtools_position_sort_strandseq_reads:
    input:
        os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            'temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam')
    output:
        temp(os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            'temp/sort/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.sam.bam')
        )
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule hac_mark_duplicate_reads_strandseq:
    input:
        rules.hac_samtools_position_sort_strandseq_reads.output[0]
    output:
        os.path.join(
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup.sam.bam')
    log:
        os.path.join('log',
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup.log')
    benchmark:
        os.path.join('run',
            'output/alignments/strandseq_to_phased_assembly',
            'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
            '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
            '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup' + '.t{}.rsrc'.format(config['num_cpu_low']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input} {output} &> {log}'


def collect_haploid_assembly_strandseq_alignments(wildcards, glob_collect=False):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    source_path = os.path.join(
        'output/alignments/strandseq_to_phased_assembly',
        'strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/{pol_reads}',
        '{hap_reads}-{hap_assembler}.{hap}.{pol_pass}',
        '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup.sam.bam{ext}'
    )

    individual, project, platform_spec = wildcards.sseq_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')

    if glob_collect:
        import glob
        source_path = source_path.replace('{ext}', '*')
        source_path = source_path.replace('{lib_id}', '*')
        wildcard_values = dict(wildcards)
        wildcard_values['individual'] = individual
        wildcard_values['project'] = project
        wildcard_values['platform'] = platform
        wildcard_values['spec'] = spec
        pattern = source_path.format(**wildcard_values)

        bam_files = glob.glob(pattern)
        if not bam_files:
            raise RuntimeError('collect_haploid_assembly_strandseq_alignments: no files collected with pattern {}'.format(pattern))

    else:
        requests_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sseq_reads).output[0]

        # this is a bit tricky given that there are different varieties of Strand-seq libraries
        glob_pattern = '_'.join([individual, project, platform + '-{spec,[0-9a-z]+}', '{lib_id}'])

        if wildcards.sseq_reads in CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES:
            glob_pattern += '_{run_id,[A-Z0-9]+}_1.request'
        else:
            glob_pattern += '_1.request'
        search_path = os.path.join(requests_dir, glob_pattern)

        checkpoint_wildcards = glob_wildcards(search_path)

        bam_files = expand(
            source_path,
            hap_assm_mode=wildcards.hap_assm_mode,
            var_caller=wildcards.var_caller,
            qual=wildcards.qual,
            gq=wildcards.gq,
            reference=wildcards.reference,
            vc_reads=wildcards.vc_reads,
            individual=individual,
            sseq_reads=wildcards.sseq_reads,
            project=project,
            platform=platform,
            spec=spec,
            lib_id=checkpoint_wildcards.lib_id,
            hap_reads=wildcards.hap_reads,
            hap_assembler=wildcards.hap_assembler,
            hap=wildcards.hap,
            pol_reads=wildcards.pol_reads,
            pol_pass=wildcards.pol_pass,
            ext=['', '.bai'])

    return sorted(bam_files)


rule hac_write_saarclust_config_file:
    """
    """
    input:
        setup_ok = rules.install_rlib_saarclust.output.check,
        reference = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.fasta',
        ref_idx = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.fasta.fai',
        strandseq_reads = 'input/fastq/{sseq_reads}.fofn',
        bam =  collect_haploid_assembly_strandseq_alignments
    output:
        cfg = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust.config',
        input_dir = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust.input'
    params:
        min_contig_size = config['min_contig_size'],
        min_region_size = config['min_region_to_order'],
        bin_size = config['bin_size'],
        step_size = config['step_size'],
        prob_threshold = config['prob_threshold'],
        init_clusters = config['init_clusters'],
        desired_clusters = config.get('desired_clusters', 0),
        min_mapq = config.get('min_mapq', 0)
    run:
        import os

        try:
            validate_checkpoint_output(input.bam)
            bam_files = input.bam
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            bam_files = collect_haploid_assembly_strandseq_alignments(wildcards, glob_collect=True)

        outfolder = os.path.dirname(bam_files[0])

        # Note to self: duplicated code somewhat needed atm
        # if non-default parameters are set for a sample.
        # Current upside: Snakemake can detect missing config values
        # in a dry run if accessed in "params" section, but not if
        # accessed in a "run" block -> fails early
        # Get rid of this via...
        # TODO: implement a (config) parameter load function to hide this
        min_contig_size = str(params.min_contig_size)
        min_region_size = str(params.min_region_size)
        bin_size = str(params.bin_size)
        step_size = str(params.step_size)
        prob_threshold = str(params.prob_threshold)
        init_clusters = str(params.init_clusters)
        desired_clusters = str(params.desired_clusters)
        min_mapq = str(params.min_mapq)

        individual = wildcards.sseq_reads.split('_')[0]
        non_default_params = config.get('sample_non_default_parameters', dict())
        if individual in non_default_params:
            sample_non_defaults = non_default_params[individual]
            use_non_defaults = True
            if 'use_only_in' in sample_non_defaults:
                try:
                    sample_non_defaults = sample_non_defaults['use_only_in']['hac_write_saarclust_config_file']
                except KeyError:
                    use_non_defaults = False
            if use_non_defaults:
                min_contig_size = str(sample_non_defaults.get('min_contig_size', min_contig_size))
                min_region_size = str(sample_non_defaults.get('min_region_size', min_region_size))
                bin_size = str(sample_non_defaults.get('bin_size', bin_size))
                step_size = str(sample_non_defaults.get('step_size', step_size))
                prob_threshold = str(sample_non_defaults.get('prob_threshold', prob_threshold))
                init_clusters = str(sample_non_defaults.get('init_clusters', init_clusters))
                desired_clusters = str(sample_non_defaults.get('desired_clusters', desired_clusters))
                min_mapq = str(sample_non_defaults.get('min_mapq', min_mapq))

        config_rows = [
            '[SaaRclust]',
            'min.contig.size = ' + min_contig_size,
            'min.region.to.order = ' + min_region_size,
            'bin.size = ' + bin_size,
            'step.size = ' + step_size,
            'prob.th = ' + prob_threshold,
            'pairedReads = TRUE',
            'store.data.obj = TRUE',
            'reuse.data.obj = FALSE',
            'num.clusters = ' + init_clusters,
            'bin.method = "dynamic"',
            'ord.method = "greedy"',
            'assembly.fasta = "' + input.reference + '"',
            'concat.fasta = FALSE',
            'remove.always.WC = TRUE',
            'mask.regions = FALSE'
        ]

        if int(config['git_commit_version']) > 7:
            config_rows.append('desired.num.clusters = ' + desired_clusters)

        if int(config['git_commit_version']) > 8:
            config_rows.append('min.mapq = ' + min_mapq)

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


checkpoint run_saarclust_haploid_assembly_clustering:
    input:
        cfg = rules.hac_write_saarclust_config_file.output.cfg,
        fofn = rules.hac_write_saarclust_config_file.output.input_dir
    output:
        dir_fasta = directory('output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust_run/clustered_assembly'),
        dir_data = directory('output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust_run/data'),
        dir_plots = directory('output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust_run/plots'),
        cfg = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust_run/SaaRclust.config',
    log:
        'log/output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust.log'
    benchmark:
        'run/output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust.rsrc'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('run_saarclust.R'),
        out_folder = lambda wildcards, output: os.path.dirname(output.cfg),
        in_folder = lambda wildcards, input: load_fofn_file(input)
    shell:
        '{params.script_exec} {input.cfg} {params.in_folder} {params.out_folder} &> {log} '



def collect_haploid_clustered_fasta_sequences(wildcards, glob_collect=False):
    """
    """
    source_path = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}/saarclust_run/clustered_assembly/{sequence}.fasta'

    if glob_collect:
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            raise RuntimeError('collect_haploid_clustered_fasta_sequences: no files collected with pattern {}'.format(pattern))

    else:

        # this output folder is the /clustered_assembly subfolder
        seq_output_dir = checkpoints.run_saarclust_haploid_assembly_clustering.get(**wildcards).output.dir_fasta
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

        fasta_files = expand(
            source_path,
            hap_assm_mode=wildcards.hap_assm_mode,
            var_caller=wildcards.var_caller,
            qual=wildcards.qual,
            gq=wildcards.gq,
            reference=wildcards.reference,
            vc_reads=wildcards.vc_reads,
            sseq_reads=wildcards.sseq_reads,
            pol_reads=wildcards.pol_reads,
            hap_reads=wildcards.hap_reads,
            hap_assembler=wildcards.hap_assembler,
            hap=wildcards.hap,
            pol_pass=wildcards.pol_pass,
            sequence=checkpoint_wildcards.sequence
        )

    return fasta_files


rule hac_write_haploid_assembly_clustered_fasta_fofn:
    """
    Local rule with minimal overhead - properly collect checkpoint output
    """
    input:
        fasta = collect_haploid_clustered_fasta_sequences
    output:
        fofn = 'output/diploid_assembly/strandseq_{hap_assm_mode}/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}/polishing/{pol_reads}/clustering/{hap_reads}-{hap_assembler}.{hap}.{pol_pass}' + '.scV{}.fofn'.format(config['git_commit_version'])
    run:
        try:
            validate_checkpoint_output(input.fasta)
            fasta_files = input.fasta
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            fasta_files = collect_haploid_clustered_fasta_sequences(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fasta_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_haploid_assembly_fasta_clusters:
    input:
        fofn = rules.hac_write_haploid_assembly_clustered_fasta_fofn.output.fofn
    output:
        rules.hac_write_haploid_assembly_clustered_fasta_fofn.output.fofn.replace('.fofn', '.fasta')
    conda:
        '../environment/conda/conda_shelltools.yml'
    params:
        fasta_clusters = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fasta_clusters} > {output}'
