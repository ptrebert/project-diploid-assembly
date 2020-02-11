
localrules: master_preprocess_input,
            write_fastq_input_parts_fofn,
            write_bam_input_parts_fofn,
            merge_strandseq_libraries


rule master_preprocess_input:
    input:
        []


def collect_fastq_input_parts(wildcards, glob_collect=False):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    import os
    subfolder = 'fastq/partial/parts'
    sample = wildcards.mrg_sample

    if glob_collect:
        import glob
        pattern = os.path.join('input', subfolder, sample + '.part*.fastq.gz')
        fastq_parts = glob.glob(pattern)

        if not fastq_parts:
            raise RuntimeError('collect_fastq_input_parts: no files collected with pattern {}'.format(pattern))

    else:

        requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

        base_path = os.path.join('input', subfolder)
        request_path = os.path.join(base_path, 'requests')

        checkpoint_wildcards = glob_wildcards(os.path.join(request_path, sample + '.{part_num}.request'))

        fastq_parts = expand(
            os.path.join(base_path, sample + '.{part_num}.fastq.gz'),
            part_num=checkpoint_wildcards.part_num
        )

    return fastq_parts


rule write_fastq_input_parts_fofn:
    input:
        req_dir = 'input/fastq/partial/parts/requests',
        fastq_parts = collect_fastq_input_parts
    output:
        fofn = 'input/fastq/complete/{mrg_sample}_1000.fofn'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    run:
        try:
            validate_checkpoint_output(input.fastq_parts)
            fastq_files = input.fastq_parts
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            fastq_files = collect_fastq_input_parts(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fastq_files):
                _ = dump.write(file_path + '\n')


rule merge_fastq_input_parts:
    input:
        fofn = 'input/fastq/complete/{mrg_sample}_1000.fofn'
    output:
        'input/fastq/complete/{mrg_sample}_1000.fastq.gz'
    log:
        'log/input/fastq/complete/{mrg_sample}_1000.merge.log'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    conda:
        '../environment/conda/conda_shelltools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.fastq_parts} > {output} 2> {log}'


def collect_pacbio_bam_input_parts(wildcards, glob_collect=False):
    """
    :param wildcards:
    :return:
    """
    import os
    subfolder = 'bam/partial/parts'
    sample = wildcards.mrg_sample

    if glob_collect:
        import glob
        pattern = os.path.join('input', subfolder, sample + '.part*.pbn.bam')
        bam_parts = glob.glob(pattern)

        if not bam_parts:
            raise RuntimeError('collect_pacbio_bam_input_parts: no files collected with pattern {}'.format(pattern))

    else:
        requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

        base_path = os.path.join('input', subfolder)
        request_path = os.path.join(base_path, 'requests')

        checkpoint_wildcards = glob_wildcards(os.path.join(request_path, sample + '.{part_num}.request'))

        bam_parts = expand(
            os.path.join(base_path, sample + '.{part_num}.pbn.bam'),
            part_num=checkpoint_wildcards.part_num
        )

    return bam_parts


rule write_bam_input_parts_fofn:
    input:
        req_dir = 'input/bam/partial/parts/requests',
        bam_parts = collect_pacbio_bam_input_parts
    output:
        fofn = 'input/bam/complete/{mrg_sample}_1000.pbn.fofn'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_PBN_INPUT_SAMPLES
    run:
        try:
            validate_checkpoint_output(input.bam_parts)
            bam_files = input.bam_parts
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            bam_files = collect_pacbio_bam_input_parts(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(bam_files):
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
        mrg_sample = CONSTRAINT_PARTS_PBN_INPUT_SAMPLES
    conda:
         '../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 if (attempt <= 1 and '-ccs' in wildcards.mrg_sample) else 24 * attempt
    params:
        bam_parts = lambda wildcards, input: load_fofn_file(input, prefix=' -in ', sep=' -in ')
    shell:
        'bamtools merge {params.bam_parts} -out {output}'


rule chs_child_filter_to_100x:
    """
    This one sample has ~200x coverage, and cannot be processed by flye
    Hard-code for now as this is not required for any other input sample
    """
    input:
        'input/bam/complete/HG00514_hgsvc_pbsq2-clr_1000.pbn.bam'
    output:
        'input/bam/complete/HG00514_hgsvc_pbsq2-clr_0526.pbn.bam'
    log:
        'log/input/bam/complete/HG00514_hgsvc_pbsq2-clr_0526.sampling.log'
    benchmark:
        'run/input/bam/complete/HG00514_hgsvc_pbsq2-clr_0526.sampling.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'bamtools filter -length ">=30000" -in {input} -out {output} &> {log}'


def collect_strandseq_libraries(wildcards, glob_collect=False):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    import os
    source_path = 'input/fastq/strand-seq/{sts_reads}/{lib_id}.fastq.gz'

    if glob_collect:
        import glob
        pattern = source_path.replace('{lib_id}', '*')
        pattern = pattern.format(**dict(wildcards))
        sseq_fastq = glob.glob(pattern)

        if not sseq_fastq:
            raise RuntimeError('collect_strandseq_libraries: no files collected with pattern {}'.format(pattern))

    else:
        checkpoint_dir = checkpoints.create_bioproject_download_requests.get(sts_reads=wildcards.sts_reads).output[0]

        glob_pattern = os.path.join(checkpoint_dir, '{lib_id}.request')

        checkpoint_wildcards = glob_wildcards(glob_pattern)

        sseq_fastq = expand(
            source_path,
            sts_reads=wildcards.sts_reads,
            lib_id=checkpoint_wildcards.lib_id
            )

    return sseq_fastq


rule merge_strandseq_libraries:
    """
    To have a simple way of incorporating the sts_reads
    wildcard into the workflow, create this file listing
    to be referred to downstream
    """
    input:
        sseq_libs = collect_strandseq_libraries
    output:
        'input/fastq/strand-seq/{sts_reads}.fofn'
    run:
        try:
            validate_checkpoint_output(input.sseq_libs)
            sseq_fastq = input.sseq_libs
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            sseq_fastq = collect_strandseq_libraries(wildcards, glob_collect=True)

        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(sseq_fastq)))
