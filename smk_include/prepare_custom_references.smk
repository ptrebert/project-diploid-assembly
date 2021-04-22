
localrules: master_prepare_custom_references,
            write_saarclust_config_file,
            write_reference_fasta_clusters_fofn,


rule master_prepare_custom_references:
    input:
        []


def collect_strandseq_merge_files(wildcards, glob_collect=False):
    """
    """

    source_path = 'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{individual}_{project}_{platform}-{spec}_{lib_id}_{run_id}.filt.sam.bam'

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
            raise RuntimeError('collect_strandseq_merge_files: no files collected with pattern {}'.format(pattern))

    else:
        requests_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sseq_reads).output[0]
        search_pattern = '_'.join([individual, project, platform + '-{spec}', lib_id, '{run_id}', '1'])

        search_path = os.path.join(requests_dir, search_pattern + '.request')

        checkpoint_wildcards = glob_wildcards(search_path)

        bam_files = expand(
            source_path,
            zip,
            reference=[wildcards.reference, wildcards.reference],
            individual=[individual, individual],
            sseq_reads=[wildcards.sseq_reads, wildcards.sseq_reads],
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
        fofn = 'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.fofn'
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
        temp('output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.rsrc'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES
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
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{library_id}.filt.sam.bam'
    output:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{library_id}.mrg.sam.bam'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES
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
        strandseq_reads = 'input/fastq/{sseq_reads}.fofn',
        bam = collect_strandseq_alignments  # from module aux_utilities
    output:
        cfg = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sseq_reads}/saarclust.config',
        input_dir = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sseq_reads}/saarclust.input'
    params:
        saarclust = lambda wildcards, input: load_saarclust_params(wildcards, input, 'squashed')
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

        with open(output.cfg, 'w') as dump:
            _ = dump.write(params.saarclust)

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


checkpoint run_saarclust_assembly_clustering:
    input:
        cfg = rules.write_saarclust_config_file.output.cfg,
        fofn = rules.write_saarclust_config_file.output.input_dir
    output:
        dir_fasta = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/clustered_assembly'),
        dir_data = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/data'),
        dir_plots = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/plots'),
        cfg = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/SaaRclust.config',
    log:
        'log/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/saarclust.log'
    benchmark:
        'run/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/saarclust.rsrc'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('run_saarclust.R'),
        out_folder = lambda wildcards, output: os.path.dirname(output.cfg),
        in_folder = lambda wildcards, input: load_fofn_file(input),
        version = config['git_commit_saarclust']
    shell:
        '{params.script_exec} {input.cfg} {params.in_folder} {params.out_folder} {params.version} &> {log} '


def collect_clustered_fasta_sequences(wildcards, glob_collect=False):
    """
    """
    import sys
    debug = bool(config.get('show_debug_messages', False))
    func_name = '\nchk::agg::collect_clustered_fasta_sequences: {}\n'

    if debug:
        sys.stderr.write(func_name.format('wildcards ' + str(wildcards)))

    source_path = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/clustered_assembly/{sequence}.fasta'

    if glob_collect:
        if debug:
            sys.stderr.write(func_name.format('called w/ glob collect'))
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            if debug:
                sys.stderr.write(func_name.format('glob collect failed'))
            raise RuntimeError('collect_clustered_fasta_sequences: no files collected with pattern {}'.format(pattern))

    else:
        if debug:
            sys.stderr.write(func_name.format('called w/ chk::get'))
        from snakemake.exceptions import IncompleteCheckpointException as ICE

        strandseq_reads = wildcards.sseq_reads
        nhr_assembly = wildcards.reference

        try:
            # this output folder is the /clustered_assembly subfolder
            seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=nhr_assembly, sseq_reads=strandseq_reads).output.dir_fasta
            checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

            fasta_files = expand(
                source_path,
                reference=nhr_assembly,
                sseq_reads=strandseq_reads,
                sequence=checkpoint_wildcards.sequence
            )
        except ICE as ice:
            if debug:
                sys.stderr.write(func_name.format('chk::get raised SMK::ICE'))
            try:
                fasta_files = collect_clustered_fasta_sequences(wildcards, glob_collect=True)
            except RuntimeError:
                if debug:
                    sys.stderr.write(func_name.format('glob collect failed - re-raising SMK::ICE'))
                raise ice
            if debug:
                sys.stderr.write('glob collect success - SMK::ICE raised in error')

    return fasta_files


rule write_reference_fasta_clusters_fofn:
    """
    Local rule with minimal overhead - properly collect checkpoint output
    """
    input:
        fasta = collect_clustered_fasta_sequences
    output:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.fofn'
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


rule check_max_cluster_size:
    """
    This is a sanity check to avoid that downstream contig alignment tasks fail
    due to squashed assembly clusters being too large to be processed (restriction of SAM/BAM)
    see:
    https://github.com/lh3/minimap2/issues/440#issuecomment-508052956
    """
    input:
        'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.fofn'
    output:
        'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.size.ok'
    run:
        max_seq_len = 268435456
        with open(input[0], 'r') as fasta_list:
            for fasta_file in fasta_list:
                seq_len = 0
                with open(fasta_file.strip(), 'r') as fasta:
                    for line in fasta:
                        if line.startswith('>'):
                            continue
                        seq_len += len(line.strip())
                if seq_len > max_seq_len:
                    raise ValueError('Squashed assembly cluster too large ({} / max {}): {}'.format(seq_len, max_seq_len, fasta_file))

        with open(output[0], 'w') as check_ok:
            pass


rule merge_reference_fasta_clusters:
    input:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{hap_reads}_nhr-{assembler}.clusters.fofn',
        size_ok = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{hap_reads}_nhr-{assembler}.clusters.size.ok'
    output:
        expand('output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{version}-{{assembler}}.fasta',
                version=config['git_commit_version'])
    conda:
        '../environment/conda/conda_shelltools.yml'
    params:
        fasta_clusters = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fasta_clusters} > {output}'
