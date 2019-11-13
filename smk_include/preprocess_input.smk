
include: 'handle_data_download.smk'
include: 'link_data_sources.smk'

localrules: master_preprocess_input


rule master_preprocess_input:
    input:


def collect_fastq_input_parts(wildcards):

    subfolder = 'fastq/partial/parts'

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

    base_path = os.path.join('input', subfolder)
    request_path = os.path.join(base_path, 'requests')

    sample = wildcards.sample

    checkpoint_wildcards = glob_wildcards(os.path.join(request_path, sample + '.{part_num}.request'))

    fastq_parts = expand(
        os.path.join(base_path, sample + '.{part_num}.fastq.gz'),
        part_num=checkpoint_wildcards.part_num
    )

    return fastq_parts


rule write_fastq_input_parts_fofn:
    input:
        collect_fastq_input_parts
    output:
        fofn = 'input/fastq/complete/{sample}_1000.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    wildcard_constraints:
        sample = '(' + '|'.join(config['partial_samples']) + ')'
    run:
        fastq_parts = collect_fastq_input_parts(wildcards)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fastq_parts):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for FASTQ merge part, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_fastq_input_parts:
    input:
        fofn = 'input/fastq/complete/{sample}_1000.fofn'
    output:
        'input/fastq/complete/{sample}_1000.fastq.gz'
    wildcard_constraints:
        sample = '(' + '|'.join(config['partial_samples']) + ')'
    log:
        'log/input/fastq/complete/{sample}_1000.merge.log'
    resources:
        runtime_hrs = 4
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fastq_parts} > {output} 2> {log}'


def collect_strandseq_libraries(wildcards):

    sample = wildcards.sample
    bioproject = config['strandseq_to_bioproject'][sample]
    individual = sample.split('_')[0]

    checkpoint_dir = checkpoints.create_bioproject_download_requests.get(individual=individual, bioproject=bioproject).output[0]

    glob_pattern = os.path.join(checkpoint_dir, '{lib_id}.request')

    checkpoint_wildcards = glob_wildcards(glob_pattern)

    checkpoint_root = os.path.split(os.path.split(checkpoint_dir)[0])[0]
    fastq_input = os.path.join(checkpoint_root, '{individual}_{bioproject}', '{lib_id}.fastq.gz')

    fastq_files = expand(
        fastq_input,
        individual=individual,
        bioproject=bioproject,
        lib_id=checkpoint_wildcards.lib_id
        )

    return fastq_files


rule merge_strandseq_libraries:
    """
    TODO
    This creates a mock file for
    the time being to keep the pipeline
    in a coherent state while the strand-seq
    steps are being implemented
    """
    input:
        collect_strandseq_libraries
    output:
        'input/fastq/complete/{sample}_sseq.fastq.gz'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    run:
        fastq_files = collect_strandseq_libraries(wildcards)
        with open(output[0], 'w') as dump:
            for fastq in sorted(fastq_files):
                dump.write(fastq + '\n')


def collect_pacbio_bam_input_parts(wildcards):

    subfolder = 'bam/partial/parts'

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

    base_path = os.path.join('input', subfolder)
    request_path = os.path.join(base_path, 'requests')

    sample = wildcards.sample

    checkpoint_wildcards = glob_wildcards(os.path.join(request_path, sample + '.{part_num}.request'))

    bam_parts = expand(
        os.path.join(base_path, sample + '.{part_num}.pbn.bam'),
        part_num=checkpoint_wildcards.part_num
    )

    return bam_parts


rule write_bam_input_parts_fofn:
    input:
        collect_pacbio_bam_input_parts
    output:
        fofn = 'input/bam/complete/{sample}_1000.pbn.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    wildcard_constraints:
        sample = '(' + '|'.join(config['partial_samples']) + ')'
    run:
        bam_parts = collect_pacbio_bam_input_parts(wildcards)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(bam_parts):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for PBN-BAM merge part, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_pacbio_native_bams:
    input:
        fofn = 'input/bam/complete/{sample}_1000.pbn.fofn'
    output:
        'input/bam/complete/{sample}_1000.pbn.bam'
    log:
        'log/input/bam/complete/{sample}_1000.mrg.log'
    benchmark:
        'run/input/bam/complete/{sample}_1000.mrg.rsrc'
    resources:
        runtime_hrs = 16
    params:
        bam_parts = lambda wildcards, input: load_fofn_file(input, prefix=' -in ', sep=' -in ')
    shell:
        'bamtools merge {params.bam_parts} -out {output}'
