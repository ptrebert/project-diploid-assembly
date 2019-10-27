
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'run_alignments.smk'

localrules: master_prepare_custom_references, \
            install_rlib_saarclust, \
            merge_mono_dinucleotide_fraction, \
            merge_reference_fasta_clusters

rule master_prepare_custom_references:
    input:


rule compute_wtdbg_squashed_assembly_layout:
    input:
        fastq = 'input/fastq/complete/{sample}_1000.fastq.gz'
    output:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
    log: 'log/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.log',
    benchmark: 'run/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 409600
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.sample],
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads}' \
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_squashed_assembly_consensus:
    input:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz'
    output:
        squashed_assembly = 'references/assemblies/{sample}_sqa-wtdbg.fasta'
    log: 'log/references/assemblies/{sample}_sqa-wtdbg.consensus.log'
    benchmark: 'run/references/assemblies/{sample}_sqa-wtdbg.consensus.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.squashed_assembly} &> {log}'


rule filter_squashed_assembly_by_size:
    input:
        'references/assemblies/{sample}_{assembly_type}-{assembler}.fasta'
    output:
        fasta = 'references/assemblies/{sample}_{assembly_type}-{assembler}-100kb.fasta',
        stats = 'output/statistics/assemblies/{sample}_{assembly_type}-{assembler}-100kb.stats.tsv'
    log:
        'log/references/assemblies/{sample}_{assembly_type}-{assembler}-100kb.log'
    wildcard_constraints:
        assembly_type = '(sqa|scV[0-9])',
    params:
        scriptdir = config['script_dir'],
        min_contig_size = config['min_contig_size'],
    shell:
        '{params.scriptdir}/filter_squashed_assembly.py --debug --input-fasta {input}' \
            ' --output-fasta {output.fasta} --output-metrics {output.stats}' \
            ' --min-size {params.min_contig_size} &> {log}'

# =================================
# ====== STRAND_SEQ PART ==========
# =================================
# Below this point: SaarClust step
# cluster squashed assembly contigs
# based on strand-seq information

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
         'output/check_files/R_setup/saarclust.ok'
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_saarclust']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_saarclust.R {params.version} &> {output}'


rule write_saarclust_config_file:
    input:
        setup_ok = 'output/check_files/R_setup/saarclust.ok',
        reference = 'references/assemblies/{reference}.fasta',
    output:
        cfg = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.config'
    params:
        min_contig_size = config['min_contig_size'],
        step_size = lambda wildcards: config['saarclust_step_size'][wildcards.reference.rsplit('-', 1)[0]],
        zlimit = lambda wildcards: config['saarclust_z_limit'][wildcards.reference.rsplit('-', 1)[0]]
    run:

        config_rows = [
            '[SaaRclust]',
            'min.contig.size = ' + str(params.min_contig_size),
            'pairedReads = TRUE',
            'bin.size = ' + str(params.min_contig_size),
            'store.data.obj = TRUE',
            'reuse.data.obj = TRUE',
            'num.clusters = 100',
            'alpha = 0.1',
            'best.prob = 1',
            'prob.th = 0',
            'ord.method = "TSP"',
            'assembly.fasta = "' + input.reference + '"',
            'concat.fasta = TRUE',
            'z.limit = ' + str(params.zlimit),
            'remove.always.WC = FALSE'
        ]

        if int(params.step_size) > 0:
            config_rows.append('step.size = ' + str(params.step_size))

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')


checkpoint run_saarclust_assembly_clustering:
    input:
        cfg = 'output/saarclust/config_files/{reference}/{sts_reads}/saarclust.config',
        bam = collect_strandseq_alignments,
    output:
        directory('output/saarclust/results/{reference}/{sts_reads}/clustered_assembly'),
        directory('output/saarclust/results/{reference}/{sts_reads}/data'),
        directory('output/saarclust/results/{reference}/{sts_reads}/plots'),
        'output/saarclust/results/{reference}/{sts_reads}/SaaRclust.config'
    log:
        'log/output/saarclust/results/{reference}/{sts_reads}/saarclust.log'
    benchmark:
        'run/output/saarclust/results/{reference}/{sts_reads}/saarclust.run'
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192
    params:
        script_dir = config['script_dir'],
        bam_folder = lambda wildcards, input: os.path.dirname(input.bam[0]),
        out_folder = lambda wildcards, output: os.path.dirname(output[3])
    shell:
        '{params.script_dir}/run_saarclust.R {input.cfg} {params.bam_folder} {params.out_folder} &> {log}'


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
