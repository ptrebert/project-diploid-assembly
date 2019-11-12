
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'run_assemblies.smk'
include: 'run_alignments.smk'

localrules: master_prepare_custom_references, \
            write_saarclust_config_file, \
            write_strandseq_merge_fofn, \
            write_reference_fasta_clusters_fofn


rule master_prepare_custom_references:
    input:


rule install_rlib_saarclust:
    output:
         'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust'])
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_saarclust']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_saarclust.R {params.version} &> {output}'


def collect_strandseq_merge_files(wildcards):
    """
    """
    individual = wildcards.individual
    bioproject = wildcards.bioproject
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    # This lookup is necessary when a (SaaR)-clustered assembly is processed,
    # because in this case the "sts_reads" wildcard is being used together
    # with the assembly name
    if bioproject.endswith('_sseq'):
        sts_reads = bioproject
        sts_lookup = sts_reads.rsplit('_', 1)[0]
        bioproject = config['strandseq_to_bioproject'][sts_lookup]

    requests_dir = checkpoints.create_bioproject_download_requests.get(individual=individual, bioproject=bioproject).output[0]

    search_pattern = '_'.join([individual, project, '{spec}', lib_id, '{run_id}', '1'])

    search_path = os.path.join(requests_dir, search_pattern + '.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/aln/{individual}_{project}_{spec}_{lib_id}_{run_id}.filt.sam.bam',
        zip,
        reference=[wildcards.reference, wildcards.reference],
        individual=[individual, individual],
        bioproject=[bioproject, bioproject],
        project=[project, project],
        spec=checkpoint_wildcards.spec,
        lib_id=[lib_id, lib_id],
        run_id=checkpoint_wildcards.run_id)

    return bam_files


rule write_strandseq_merge_fofn:
    input:
        bams = collect_strandseq_merge_files
    output:
        fofn = 'output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.fofn'
    run:
        bam_files = collect_strandseq_merge_files(wildcards)
        if len(bam_files) != 2:
            raise AssertionError('Missing merge partner for strand-seq BAM files: {}'.format(bam_files))

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
        temp('output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/mrg/{individual}_{project}_{platform}-npe_{lib_id}.mrg.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/mrg{individual}_{project}_{platform}-npe_{lib_id}.mrg.rsrc'
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
    threads: config['num_cpu_low']
    resources:
            mem_per_cpu_mb = 512,
            mem_total_mb = config['num_cpu_low'] * 512
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
        reference = 'output/reference_assembly/squashed/{reference}.fasta',
        ref_idx = 'output/reference_assembly/squashed/{reference}.fasta.fai',
        strandseq_reads = 'input/fastq/complete/{sts_reads}.fastq.gz',
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

# === DEBUG ====
# SaaRclust is currently detected as failing by Snakemake for unclear
# reasons; since the behavior cannot be reproduced on the command line,
# manually generated SaaRclust output has to be copied in the appropriate
# place
# To protect this manually generated output, the SaaRclust checkpoint produces
# protected output, and the subsequent merge operations use ancient()

checkpoint run_saarclust_assembly_clustering:
    input:
        cfg = ancient(rules.write_saarclust_config_file.output.cfg),
        fofn = ancient(rules.write_saarclust_config_file.output.input_dir)
    output:
        dir_fasta = protected(directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/clustered_assembly')),
        dir_data = protected(directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/data')),
        dir_plots = protected(directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/plots')),
        cfg = protected('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/SaaRclust.config'),
    log:
        'log/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/saarclust.log'
    benchmark:
        'run/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sts_reads}/saarclust.run'
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192
    params:
        script_dir = config['script_dir'],
        out_folder = lambda wildcards, output: os.path.dirname(output.cfg),
        in_folder = lambda wildcards, input: load_fofn_file(input)
    shell:
        '{params.script_dir}/run_saarclust.R {input.cfg} {params.in_folder} {params.out_folder} &> {log} ' \
            ' && echo "SaaRclust success" &>> {log} ' \
            ' || echo "SaaRclust fail" &>> {log}'


def collect_clustered_fasta_sequences(wildcards):
    """
    """
    strandseq_reads = wildcards.sts_reads
    sqa_assembly = wildcards.reference

    # this output folder is the /clustered_assembly subfolder
    seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=sqa_assembly, sts_reads=strandseq_reads).output.dir_fasta
    checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

    cluster_fasta = expand(
        os.path.join(seq_output_dir, '{sequence}.fasta'),
        reference=sqa_assembly,
        sts_reads=strandseq_reads,
        sequence=checkpoint_wildcards.sequence
    )

    return sorted(cluster_fasta)


rule write_reference_fasta_clusters_fofn:
    """
    Local rule with minimal overhead - properly collect checkpoint output
    """
    input:
        fasta = ancient(collect_clustered_fasta_sequences)
    output:
       fofn = 'output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{reference}.clusters.fofn'
    run:
        # following example as above for merge strand-seq BAM files
        fasta_files = collect_clustered_fasta_sequences(wildcards)

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
        fofn = ancient('output/reference_assembly/clustered/temp/saarclust/{sts_reads}/{sample}_sqa-{assembler}.clusters.fofn')
    output:
        expand('output/reference_assembly/clustered/{{sts_reads}}/{{sample}}_scV{version}-{{assembler}}.fasta',
                version=config['git_commit_version'])
    params:
        fasta_clusters = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fasta_clusters} > {output}'
