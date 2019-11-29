
localrules: master_preprocess_input


rule master_preprocess_input:
    input:


def collect_fastq_input_parts(wildcards):

    subfolder = 'fastq/partial/parts'

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

    base_path = os.path.join('input', subfolder)
    request_path = os.path.join(base_path, 'requests')

    sample = wildcards.mrg_sample

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
        fofn = 'input/fastq/complete/{mrg_sample}_1000.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_fastq_samples']) + ')'
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
        fofn = 'input/fastq/complete/{mrg_sample}_1000.fofn'
    output:
        'input/fastq/complete/{mrg_sample}_1000.fastq.gz'
    log:
        'log/input/fastq/complete/{mrg_sample}_1000.merge.log'
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_fastq_samples']) + ')'
    resources:
        runtime_hrs = 4
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fastq_parts} > {output} 2> {log}'


def collect_pacbio_bam_input_parts(wildcards):

    subfolder = 'bam/partial/parts'

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

    base_path = os.path.join('input', subfolder)
    request_path = os.path.join(base_path, 'requests')

    sample = wildcards.mrg_sample

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
        fofn = 'input/bam/complete/{mrg_sample}_1000.pbn.fofn'
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_pbn_samples']) + ')'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
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
        fofn = 'input/bam/complete/{mrg_sample}_1000.pbn.fofn'
    output:
        'input/bam/complete/{mrg_sample}_1000.pbn.bam'
    log:
        'log/input/bam/complete/{mrg_sample}_1000.mrg.log'
    benchmark:
        'run/input/bam/complete/{mrg_sample}_1000.mrg.rsrc'
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_pbn_samples']) + ')'
    resources:
        runtime_hrs = 16
    params:
        bam_parts = lambda wildcards, input: load_fofn_file(input, prefix=' -in ', sep=' -in ')
    shell:
        'bamtools merge {params.bam_parts} -out {output}'



def collect_strandseq_libraries(wildcards):

    checkpoint_dir = checkpoints.create_bioproject_download_requests.get(sts_reads=wildcards.sts_reads).output[0]

    glob_pattern = os.path.join(checkpoint_dir, '{lib_id}.request')

    checkpoint_wildcards = glob_wildcards(glob_pattern)

    checkpoint_root = os.path.split(os.path.split(checkpoint_dir)[0])[0]
    fastq_input = os.path.join(checkpoint_root, '{sts_reads}', '{lib_id}.fastq.gz')

    fastq_files = expand(
        fastq_input,
        sts_reads=wildcards.sts_reads,
        lib_id=checkpoint_wildcards.lib_id
        )

    return fastq_files


rule merge_strandseq_libraries:
    """
    To have a simple way of incorporating the sts_reads
    wildcard into the workflow, create this file listing
    to be referred to downstream
    """
    input:
        collect_strandseq_libraries
    output:
        'input/fastq/strand-seq/{sts_reads}.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    run:
        fastq_files = collect_strandseq_libraries(wildcards)
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(fastq_files)))
