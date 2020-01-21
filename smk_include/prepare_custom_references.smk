
localrules: master_prepare_custom_references


rule master_prepare_custom_references:
    input:
        []


def collect_strandseq_merge_files(wildcards):
    """
    """
    individual = wildcards.individual
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    assert individual in wildcards.reference and individual in wildcards.sts_reads, \
        'Wrong reference / sts_reads match: {} / {}'.format(wildcards.reference, wildcards.sts_reads)

    requests_dir = checkpoints.create_bioproject_download_requests.get(sts_reads=wildcards.sts_reads).output[0]

    search_pattern = '_'.join([individual, project, '{spec}', lib_id, '{run_id}', '1'])

    search_path = os.path.join(requests_dir, search_pattern + '.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/aln/{individual}_{project}_{spec}_{lib_id}_{run_id}.filt.sam.bam',
        zip,
        reference=[wildcards.reference, wildcards.reference],
        individual=[individual, individual],
        sts_reads=[wildcards.sts_reads, wildcards.sts_reads],
        project=[project, project],
        spec=checkpoint_wildcards.spec,
        lib_id=[lib_id, lib_id],
        run_id=checkpoint_wildcards.run_id)

    return bam_files


rule write_strandseq_merge_fofn:
    input:
        bams = collect_strandseq_merge_files
    output:
        fofn = 'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.fofn'
    run:
        import os
        bam_files = collect_strandseq_merge_files(wildcards)
        if len(bam_files) != 2:
            raise RuntimeError('Missing merge partner for strand-seq BAM files: '
                               '{}'.format(output.fofn))

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(bam_files):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for strand-seq BAM merge, but received directory: {}'.format(file_path))
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
        temp('output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.mrg.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/mrg{individual}_{project}_{platform}-npe_{lib_id}.mrg.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    params:
        merge_files = lambda wildcards, input: load_fofn_file(input)
    shell:
        'samtools merge -@ {threads} -O BAM {output} {params.merge_files} &> {log}'


rule samtools_position_sort_strandseq_reads:
    """
    Since Strand-seq alignments are small, make dedicated
    samtools sort rule with lower resource requirements
    """
    input:
        '{folder_path}/temp/mrg/{sample}.mrg.sam.bam'
    output:
        temp('{folder_path}/temp/sort/{sample}.mrg.psort.sam.bam')
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 512,
        mem_total_mb = config['num_cpu_low'] * 512
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule mark_duplicate_reads_strandseq:
    input:
        rules.samtools_position_sort_strandseq_reads.output[0]
    output:
        '{folder_path}/{sample}.mrg.psort.mdup.sam.bam'
    log:
        'log/{folder_path}/{sample}.mrg.psort.mdup.log'
    benchmark:
        'run/{folder_path}/{sample}.mrg.psort.mdup.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 512,
        mem_total_mb = config['num_cpu_low'] * 512,
        runtime_hrs = lambda wildcards, attempt: attempt + (attempt -1 ) * 2
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
        setup_ok = 'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']),
        reference = 'output/reference_assembly/non-hap-res/{reference}.fasta',
        ref_idx = 'output/reference_assembly/non-hap-res/{reference}.fasta.fai',
        strandseq_reads = 'input/fastq/strand-seq/{sts_reads}.fofn',
        bam = collect_strandseq_alignments  # from module aux_utilities
    output:
        cfg = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sts_reads}/saarclust.config',
        input_dir = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sts_reads}/saarclust.input'
    params:
        min_contig_size = config['min_contig_size'],
        bin_size = config['bin_size'],
        step_size = config['step_size'],
        prob_threshold = config['prob_threshold']
    run:
        import os
        # following same example as merge strand-seq BAMs above
        bam_files = collect_strandseq_alignments(wildcards)
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
            'num.clusters = 100',
            'bin.method = "dynamic"',
            'assembly.fasta = "' + input.reference + '"',
            'concat.fasta = TRUE',
            'remove.always.WC = TRUE',
            'mask.regions = FALSE'
        ]

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
        'run/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/saarclust.run'
    conda:
        '../environment/conda/conda_rtools.yml'
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192,
        runtime_hrs = 23
    params:
        script_dir = config['script_dir'],
        out_folder = lambda wildcards, output: os.path.dirname(output.cfg),
        in_folder = lambda wildcards, input: load_fofn_file(input)
    shell:
        '{params.script_dir}/run_saarclust.R {input.cfg} {params.in_folder} {params.out_folder} &> {log} '


def collect_clustered_fasta_sequences(wildcards):
    """
    """
    strandseq_reads = wildcards.sts_reads
    nhr_assembly = wildcards.reference

    # this output folder is the /clustered_assembly subfolder
    seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=nhr_assembly, sts_reads=strandseq_reads).output.dir_fasta
    checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

    cluster_fasta = expand(
        os.path.join(seq_output_dir, '{sequence}.fasta'),
        reference=nhr_assembly,
        sts_reads=strandseq_reads,
        sequence=checkpoint_wildcards.sequence
    )

    return sorted(cluster_fasta)


rule write_reference_fasta_clusters_fofn:
    """
    Local rule with minimal overhead - properly collect checkpoint output
    """
    input:
        fasta = collect_clustered_fasta_sequences
    output:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{reference}.clusters.fofn'
    run:
        import os
        # following example as above for merge strand-seq BAM files
        fasta_files = collect_clustered_fasta_sequences(wildcards)
        if len(fasta_files) == 0:
            raise RuntimeError('No FASTA files to merge after SaaRclust run. '
                               'SaaRclust most likely failed for sample: {}'.format(output.fofn))

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fasta_files):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for strand-seq BAM merge, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_reference_fasta_clusters:
    input:
        fofn = 'output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{sample}_nhr-{assembler}.clusters.fofn'
    output:
        expand('output/reference_assembly/clustered/{{sts_reads}}/{{sample}}_scV{version}-{{assembler}}.fasta',
                version=config['git_commit_version'])
    conda:
        '../environment/conda/conda_shelltools.yml'
    params:
        fasta_clusters = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fasta_clusters} > {output}'


# Below: create a diagnostic plot (SaaRclust) package, requires contig to reference alignment

rule dump_contig_to_reference_alignment_to_bed:
        input:
            'output/alignments/contigs_to_reference/{folder_path}/{reference}_map-to_{aln_reference}.psort.sam.bam'
        output:
            'output/alignments/contigs_to_reference/{folder_path}/{reference}_map-to_{aln_reference}.bed'
        log:
            'log/output/alignments/contigs_to_reference/{folder_path}/{reference}_map-to_{aln_reference}.dump-bed.log'
        benchmark:
            'run/output/alignments/contigs_to_reference/{folder_path}/{reference}_map-to_{aln_reference}.dump-bed.rsrc'
        conda:
            '../environment/conda/conda_biotools.yml'
        resources:
            runtime_hrs = lambda wildcards, attempt: attempt,
            mem_total_mb = 2048,
            mem_per_cpu_mb = 2048
        shell:
            'bedtools bamtobed -i {input} > {output} 2> {log}'


rule plot_saarclust_diagnostic_output:
    input:
        setup_ok = 'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']),
        ctg_ref_aln = 'output/alignments/contigs_to_reference/{folder_path}/{reference}_map-to_{aln_reference}.bed'
    output:
        clustering = 'output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.clustering.pdf',
        ordering = 'output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.ordering.pdf',
        orienting = 'output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.orienting.pdf',
    log:
       'log/output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.saarclust-diagnostics.log'
    benchmark:
        'run/output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.saarclust-diagnostics.rsrc'
    conda:
        '../environment/conda/conda_rtools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        run_script = os.path.join(config['script_dir'], 'plot_saarclust_diagnostics.R'),
        out_prefix = lambda wildcards: os.path.join(
            'output', 'plotting', 'saarclust_diagnostics', 'reference_assembly', 'clustered',
            wildcards.folder_path, wildcards.reference + '_map-to_' + wildcards.aln_reference)
    shell:
         '{params.run_script} {input.ctg_ref_aln} hg38 {params.out_prefix} &> {log}'