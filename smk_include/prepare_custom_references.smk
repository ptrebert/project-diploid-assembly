
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'run_assemblies.smk'
include: 'run_alignments.smk'

localrules: master_prepare_custom_references, \
            install_rlib_saarclust, \
            write_saarclust_config_file, \
            merge_mono_dinucleotide_fraction, \
            merge_reference_fasta_clusters

rule master_prepare_custom_references:
    input:


def collect_strandseq_merge_files(wildcards):
    """
    """
    individual = wildcards.individual
    bioproject = wildcards.bioproject
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    requests_dir = checkpoints.create_bioproject_download_requests.get(individual=individual, bioproject=bioproject).output[0]

    search_pattern = '_'.join([individual, project, '{spec}', lib_id, '{run_id}', '1'])

    search_path = os.path.join(requests_dir, search_pattern + '.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}.filt.sam.bam',
        zip,
        reference=[wildcards.reference, wildcards.reference],
        individual=[individual, individual],
        bioproject=[bioproject, bioproject],
        project=[project, project],
        spec=checkpoint_wildcards.spec,
        lib_id=[lib_id, lib_id],
        run_id=checkpoint_wildcards.run_id)

    assert len(bam_files) == 2, 'Missing merge partner: {}'.format(bam_files)

    return sorted(bam_files)


rule merge_mono_dinucleotide_fraction:
    """
    This rule is likely quite specific
    for the strand-seq data used in this
    pipeline - point of 'failure' for
    different input data
    """
    input:
        nuc_files = collect_strandseq_merge_files
    output:
        temp('output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.rsrc'
    wildcard_constraints:
        lib_id = '[A-Z0-9]+'
    threads: config['num_cpu_local']
    shell:
        'samtools merge -@ {threads} -O BAM {output} {input.nuc_files} &> {log}'


rule samtools_position_sort_strandseq_reads:
    input:
        'output/alignments/strandseq_to_reference/{subfolder}/{sample}.mrg.sam.bam'
    output:
        temp('output/alignments/strandseq_to_reference/{subfolder}/{sample}.mrg.psort.sam.bam')
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = config['num_cpu_low'] * 2048
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule mark_duplicate_reads_strandseq:
    input:
        'output/alignments/strandseq_to_reference/{subfolder}/{sample}.mrg.psort.sam.bam'
    output:
        'output/alignments/strandseq_to_reference/{subfolder}/final/{sample}.mrg.psort.mdup.sam.bam'
    log:
        'log/output/alignments/strandseq_to_reference/{subfolder}/{sample}.mrg.psort.mdup.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{subfolder}/{sample}.mrg.psort.mdup.rsrc'
    threads: config['num_cpu_low']
    shell:
        'sambamba markdup -t {threads} {input} {output} &> {log}'


rule install_rlib_saarclust:
    output:
         'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust'])
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_saarclust']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_saarclust.R {params.version} &> {output}'


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
        reference = 'references/assemblies/{reference}.fasta',
        bam = collect_strandseq_alignments
    output:
        cfg = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.config',
        input_dir = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.input',
    params:
        min_contig_size = config['min_contig_size'],
        bin_size = config['bin_size'],
        step_size = config['step_size'],
        prob_threshold = config['prob_threshold']
    run:

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
            'bin.method = "fixed"',
            'assembly.fasta = "' + input.reference + '"',
            'concat.fasta = FALSE',
            'remove.always.WC = TRUE',
            'mask.regions = FALSE'
        ]

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        outfolder = os.path.dirname(input.bam[0])
        assert os.path.isdir(outfolder), 'Invalid output folder / strand-seq alignments: {}'.format(outfolder)
        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


checkpoint run_saarclust_assembly_clustering:
    input:
        cfg = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.config',
        dir = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.input',
    output:
        dir_fasta = directory('output/saarclust/results/{reference}/{sts_reads}/clustered_assembly'),
        dir_data = directory('output/saarclust/results/{reference}/{sts_reads}/data'),
        dir_plots = directory('output/saarclust/results/{reference}/{sts_reads}/plots'),
        cfg = 'output/saarclust/results/{reference}/{sts_reads}/SaaRclust.config',
    log:
        'log/output/saarclust/results/{reference}/{sts_reads}/saarclust.log'
    benchmark:
        'run/output/saarclust/results/{reference}/{sts_reads}/saarclust.run'
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192
    params:
        script_dir = config['script_dir'],
        out_folder = lambda wildcards, output: os.path.dirname(output[3])
    run:
        with open(input.dir) as info:
            input_folder = info.readline().strip()

        exec = '{params.script_dir}/run_saarclust.R {input.cfg}'
        exec += ' {} '.format(input_folder)
        exec += '{params.out_folder} &> {log}'
        shell(exec)


def collect_clustered_fasta_sequences(wildcards):
    """
    """
    strandseq_reads = config['strandseq_to_assembly'][wildcards.reference]

    sqa_assembly = wildcards.reference + '_sqa-' + wildcards.assembler

    # this output folder is the /clustered_assembly subfolder
    seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=sqa_assembly, sts_reads=strandseq_reads).output[0]
    checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.seq'))

    cluster_fasta = expand(
        'output/saarclust/results/{reference}/{sts_reads}/clustered_assembly/{sequence}.fasta',
        reference=sqa_assembly,
        sts_reads=strandseq_reads,
        sequence=checkpoint_wildcards.sequence
        )

    return sorted(cluster_fasta)


rule merge_reference_fasta_clusters:
    input:
        fasta = collect_clustered_fasta_sequences
    output:
        expand('references/assemblies/{{reference}}_scV{version}-{{assembler}}.fasta',
                version=config['git_commit_version'])
    threads: config['num_cpu_local']
    shell:
        'cat {input} > {output}'
