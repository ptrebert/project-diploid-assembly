
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
    log:
        'input/fastq/complete/{mrg_sample}_1000.fofn.log'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_fastq_samples']) + ')'
    run:
        fastq_parts = collect_fastq_input_parts(wildcards)

        with open(log[0], 'w') as logfile:
            _ = logfile.write('{} FASTQ parts for merging\n'.format(len(fastq_parts)))

            with open(output.fofn, 'w') as dump:
                for file_path in sorted(fastq_parts):
                    if not os.path.isfile(file_path):
                        _ = logfile.write('Warning: path is not (yet) a file on the file system: {}\n'.format(file_path))
                        if os.path.isdir(file_path):
                            _ = logfile.write('Error: assumed file path is TRUE for directory: {}\n'.format(file_path))
                            raise ValueError('Expected file path for PBN-BAM merge part, but received directory: {}'.format(file_path))
                    _ = dump.write(file_path + '\n')
                    records_written += 1
            if records_written == 0:
                _ = logfile.write('Error: 0 records written to fofn file: {}\n'.format(output.fofn))
                os.unlink(output.fofn)
            else:
                _ = logfile.write('Completed: {} records written to fofn file: {}'.format(records_written, output.fofn))


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
    log:
        'input/bam/complete/{mrg_sample}_1000.pbn.fofn.log'
    wildcard_constraints:
        mrg_sample = '(' + '|'.join(config['partial_pbn_samples']) + ')'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    run:
        bam_parts = collect_pacbio_bam_input_parts(wildcards)
        records_written = 0

        with open(log[0], 'w') as logfile:
            _ = logfile.write('{} pbn.bam parts for merging\n'.format(len(bam_parts)))

            with open(output.fofn, 'w') as dump:
                for file_path in sorted(bam_parts):
                    if not os.path.isfile(file_path):
                        _ = logfile.write('Warning: path is not (yet) a file on the file system: {}\n'.format(file_path))
                        if os.path.isdir(file_path):
                            _ = logfile.write('Error: assumed file path is TRUE for directory: {}\n'.format(file_path))
                            raise ValueError('Expected file path for PBN-BAM merge part, but received directory: {}'.format(file_path))
                    _ = dump.write(file_path + '\n')
                    records_written += 1
            if records_written == 0:
                _ = logfile.write('Error: 0 records written to fofn file: {}\n'.format(output.fofn))
                os.unlink(output.fofn)
            else:
                _ = logfile.write('Completed: {} records written to fofn file: {}'.format(records_written, output.fofn))


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
        runtime_hrs = 48
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
