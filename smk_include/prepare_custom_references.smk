
localrules: master_prepare_custom_references,
            write_saarclust_config_file,
            write_reference_fasta_clusters_fofn,
            link_strandseq_monofraction_samples


rule master_prepare_custom_references:
    input:
        []


def collect_strandseq_merge_files(wildcards, glob_collect=False):
    """
    """

    source_path = 'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/aln/{individual}_{project}_{platform}-{spec}_{lib_id}_{run_id}.filt.sam.bam'

    individual = wildcards.individual
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    assert individual in wildcards.reference and individual in wildcards.sts_reads, \
        'Wrong reference / sts_reads match: {} / {}'.format(wildcards.reference, wildcards.sts_reads)

    if glob_collect:
        import glob
        pattern = source_path.replace('{run_id}', '*')
        pattern = pattern.replace('{spec}', '*')
        pattern = pattern.format(**dict(wildcards))
        bam_files = glob.glob(pattern)

        if not bam_files:
            raise RuntimeError('collect_strandseq_merge_files: no files collected with pattern {}'.format(pattern))

    else:
        requests_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sts_reads).output[0]
        search_pattern = '_'.join([individual, project, platform + '-{spec}', lib_id, '{run_id}', '1'])

        search_path = os.path.join(requests_dir, search_pattern + '.request')

        checkpoint_wildcards = glob_wildcards(search_path)

        bam_files = expand(
            source_path,
            zip,
            reference=[wildcards.reference, wildcards.reference],
            individual=[individual, individual],
            sts_reads=[wildcards.sts_reads, wildcards.sts_reads],
            project=[project, project],
            platform=[platform, platform],
            spec=checkpoint_wildcards.spec,
            lib_id=[lib_id, lib_id],
            run_id=checkpoint_wildcards.run_id)

    return bam_files


rule write_strandseq_merge_fofn:
    """
    2020-02-05
    Since this has to be executed once per merge pair,
    I don't want to make this a local rule executed on the submit node
    (presumably fixing github issue #216). However, the validate_checkpoint_output
    function indeed fails; apparently, issue #55 is still unfixed in 5.10.0
    So, implement the checkpoint functionality manually for this one... awesome
    """
    input:
        bams = collect_strandseq_merge_files
    output:
        fofn = 'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.fofn'
    wildcard_constraints:
        sts_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES,
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
            bam_files = collect_strandseq_merge_files(wildcards, glob_collect=True)

        if len(bam_files) != 2:
            raise RuntimeError('Missing merge partner for strand-seq BAM files {} / {}: '
                               '{}'.format(output.fofn, pattern, bam_files))

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(bam_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_mono_dinucleotide_fraction:
    """
    This rule is likely quite specific
    for the strand-seq data used in this
    pipeline - point of 'failure' for
    different input data
    """
    input:
        fofn = rules.write_strandseq_merge_fofn.output.fofn
    output:
        temp('output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.rsrc'
    wildcard_constraints:
        sts_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    params:
        merge_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'samtools merge -@ {threads} -O BAM {output} {params.merge_files} &> {log}'


rule link_strandseq_monofraction_samples:
    """
    Switch to copying here because sym linking seems to cause
    trouble (not recognized) on certain file systems
    """
    input:
        'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/aln/{library_id}.filt.sam.bam'
    output:
        'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{library_id}.mrg.sam.bam'
    wildcard_constraints:
        sts_reads = CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES
    shell:
        'cp {input} {output}'


rule samtools_position_sort_strandseq_reads:
    """
    Since Strand-seq alignments are small, make dedicated
    samtools sort rule with lower resource requirements
    """
    input:
        '{folder_path}/temp/mrg/{sts_library}.mrg.sam.bam'
    output:
        temp('{folder_path}/temp/sort/{sts_library}.mrg.psort.sam.bam')
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sts_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule mark_duplicate_reads_strandseq:
    input:
        rules.samtools_position_sort_strandseq_reads.output[0]
    output:
        '{folder_path}/{sts_library}.mrg.psort.mdup.sam.bam'
    log:
        'log/{folder_path}/{sts_library}.mrg.psort.mdup.log'
    benchmark:
        'run/{folder_path}/{sts_library}.mrg.psort.mdup.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sts_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input} {output} &> {log}'


rule write_saarclust_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for StrandPhaseR
    """
    input:
        setup_ok = rules.install_rlib_saarclust.output.check,
        reference = 'output/reference_assembly/non-hap-res/{reference}.fasta',
        ref_idx = 'output/reference_assembly/non-hap-res/{reference}.fasta.fai',
        strandseq_reads = 'input/fastq/{sts_reads}.fofn',
        bam = collect_strandseq_alignments  # from module aux_utilities
    output:
        cfg = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sts_reads}/saarclust.config',
        input_dir = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sts_reads}/saarclust.input'
    params:
        min_contig_size = config['min_contig_size'],
        bin_size = config['bin_size'],
        step_size = config['step_size'],
        prob_threshold = config['prob_threshold'],
        init_clusters = config['init_clusters'],
        desired_clusters = config.get('desired_clusters', None)
    run:
        import os

        try:
            validate_checkpoint_output(input.bam)
            bam_files = input.bam
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            bam_files = collect_strandseq_alignments(wildcards, glob_collect=True)

        outfolder = os.path.dirname(bam_files[0])

        config_rows = [
            '[SaaRclust]',
            'min.contig.size = ' + str(params.min_contig_size),
            'bin.size = ' + str(params.bin_size),
            'step.size = ' + str(params.step_size),
            'prob.th = ' + str(params.prob_threshold),
            'pairedReads = TRUE',
            'store.data.obj = TRUE',
            'reuse.data.obj = TRUE',
            'num.clusters = ' + str(params.init_clusters),
            'bin.method = "dynamic"',
            'assembly.fasta = "' + input.reference + '"',
            'concat.fasta = TRUE',
            'remove.always.WC = TRUE',
            'mask.regions = FALSE'
        ]

        if int(config['git_commit_version']) > 7:
            config_rows.append('desired.num.clusters = ' + str(params.desired_clusters))

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


checkpoint run_saarclust_assembly_clustering:
    input:
        cfg = rules.write_saarclust_config_file.output.cfg,
        fofn = rules.write_saarclust_config_file.output.input_dir
    output:
        dir_fasta = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/clustered_assembly'),
        dir_data = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/data'),
        dir_plots = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/plots'),
        cfg = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/SaaRclust.config',
    log:
        'log/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/saarclust.log'
    benchmark:
        'run/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/saarclust.rsrc'
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


def collect_clustered_fasta_sequences(wildcards, glob_collect=False):
    """
    """
    source_path = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/clustered_assembly/{sequence}.fasta'

    if glob_collect:
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            raise RuntimeError('collect_clustered_fasta_sequences: no files collected with pattern {}'.format(pattern))

    else:
        strandseq_reads = wildcards.sts_reads
        nhr_assembly = wildcards.reference

        # this output folder is the /clustered_assembly subfolder
        seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=nhr_assembly, sts_reads=strandseq_reads).output.dir_fasta
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

        fasta_files = expand(
            source_path,
            reference=nhr_assembly,
            sts_reads=strandseq_reads,
            sequence=checkpoint_wildcards.sequence
        )

    return fasta_files


rule write_reference_fasta_clusters_fofn:
    """
    Local rule with minimal overhead - properly collect checkpoint output
    """
    input:
        fasta = collect_clustered_fasta_sequences
    output:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{reference}.clusters.fofn'
    run:
        try:
            validate_checkpoint_output(input.fasta)
            fasta_files = input.fasta
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            fasta_files = collect_clustered_fasta_sequences(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fasta_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_reference_fasta_clusters:
    input:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{hap_reads}_nhr-{assembler}.clusters.fofn'
    output:
        expand('output/reference_assembly/clustered/{{sts_reads}}/{{hap_reads}}_scV{version}-{{assembler}}.fasta',
                version=config['git_commit_version'])
    conda:
        '../environment/conda/conda_shelltools.yml'
    params:
        fasta_clusters = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fasta_clusters} > {output}'
