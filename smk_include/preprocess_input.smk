
localrules: master_preprocess_input,
            write_fastq_input_parts_fofn,
            write_bam_input_parts_fofn,
            merge_strandseq_libraries,
            write_short_read_input_fofn,
            relink_complete_short_read_input_samples


rule master_preprocess_input:
    input:
        []


rule download_hifi_adapter_db:
    """
    This is not part of the general env setup b/c it only applies to HiFi reads
    and is likely to be deprecated soon (when PacBio has fixed the issue)
    """
    output:
        'references/annotation/HiFiAdapterFilt-master/DB/pacbio_vectors_db.nhr'
    conda:
        '../environment/conda/conda_shelltools.yml'
    shell:
        'wget -O references/downloads/HiFiAdapterFilt.zip '
            '-q https://github.com/sheinasim/HiFiAdapterFilt/archive/master.zip && '
            'unzip -q references/downloads/HiFiAdapterFilt.zip -d references/annotation/ '


rule blast_hifi_reads_adapter_contamination:
    """
    Scanning for adapter contamination follows this implementation:
        github.com/sheinasim/HiFiAdapterFilt
    BLAST output format #6:
        qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        [http://www.metagenomics.wiki/tools/blast/blastn-output-format-6]
    """
    input:
        fastq_part = 'input/fastq/{readset}/{readset}.{partnum}.fastq.gz',
        adapter_db = 'references/annotation/HiFiAdapterFilt-master/DB/pacbio_vectors_db.nhr',
    output:
        'input/fastq/{readset}/{readset}.{partnum}.blast.out',
    benchmark:
        'rsrc/input/fastq/{readset}/{readset}.{partnum}.blast.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 2 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((4096 * attempt) / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
    params:
        adapter_db_prefix = 'references/annotation/HiFiAdapterFilt-master/DB/pacbio_vectors_db',
        blast_threads = int(config['num_cpu_low']) - 1
    shell:
        'seqtk seq -C -A {input.fastq_part} | '
            'blastn -db {params.adapter_db_prefix} -query /dev/stdin -num_threads {params.blast_threads} -task blastn '
            '-reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue .01 -searchsp 1750000000000 '
            '-outfmt 6 > {output}'


rule build_hifi_read_blocklist:
    """
    Filters apply to columns "% identity" and "match length"
    """
    input:
        'input/fastq/{readset}/{readset}.{partnum}.blast.out',
    output:
        'input/fastq/{readset}/{readset}.{partnum}.contaminated-reads.txt',
    conda:
        '../environment/conda/conda_shelltools.yml'
    shell:
        "grep 'NGB0097' {input} | awk -v OFS='\t' '{{if (($2 ~ /NGB00972/ && $3 >= 97 && $4 >= 44) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}}' | sort -u > {output}"


rule select_clean_hifi_reads:
    input:
        fastq_part = 'input/fastq/{readset}/{readset}.{partnum}.fastq.gz',
        blocklist = 'input/fastq/{readset}/{readset}.{partnum}.contaminated-reads.txt',
    output:
        'input/fastq/{readset}/{readset}.{partnum}.clean-reads.fastq.gz'
    benchmark:
        'rsrc/input/fastq/{readset}/{readset}.{partnum}.clean-reads.filter.rsrc'
    wildcard_constraints:
        readset = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    resources:
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    conda:
        '../environment/conda/conda_biotools.yml'
    shell:
        'seqtk seq -C {input.fastq_part} | paste - - - - | grep -v -F -f {input.blocklist} | tr "\t" "\n" | gzip > {output}'


def collect_fastq_input_parts(wildcards, glob_collect=True, caller='snakemake'):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    import os
    import sys
    sample = wildcards.mrg_sample
    subfolder = os.path.join('fastq', sample)

    infix = ''
    try:
        check_adapter_contamination = bool(config['check_hifi_adapter_contamination'])
    except KeyError:
        sys.stderr.write('\nWarning: no key "check_hifi_adapter_contamination" in parameter config - skipping\n')
        check_adapter_contamination = False

    if 'ccs' in sample and check_adapter_contamination:
        infix = 'clean-reads.'

    if glob_collect:
        import glob
        pattern = os.path.join('input', subfolder, sample + '.part*.{}fastq.gz'.format(infix))
        fastq_parts = glob.glob(pattern)

        if not fastq_parts:
            if caller == 'snakemake':
                fastq_parts = []
                path = pattern.replace('*', '{}')
                # easy case: get number of input parts
                num_parts = count_number_of_input_parts(wildcards.mrg_sample)
                fastq_parts = [path.format(i) for i in range(1, num_parts + 1)]
            else:
                raise RuntimeError('collect_fastq_input_parts: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        request_path = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=sample).output[0]

        base_path = os.path.join('input', subfolder)

        checkpoint_wildcards = glob_wildcards(os.path.join(request_path, sample + '.{part_num}.request'))

        fastq_parts = expand(
            os.path.join(base_path, sample + '.{{part_num}}.{}fastq.gz'.format(infix)),
            part_num=checkpoint_wildcards.part_num
        )

    assert fastq_parts, 'collect_fastq_input_parts >> returned empty output: {}'.format(wildcards)
    return fastq_parts


rule write_fastq_input_parts_fofn:
    input:
        req_dir = 'input/fastq/{mrg_sample}/requests',
        fastq_parts = collect_fastq_input_parts
    output:
        fofn = 'input/fastq/{mrg_sample}_1000.fofn'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    run:
        with open(output.fofn, 'w') as dump:
            for file_path in sorted(input.fastq_parts):
                _ = dump.write(file_path + '\n')


rule merge_fastq_input_parts:
    """
    Why this compression overhead?
    pysam has issues iterating through gzip files
    that are the result of a simple concat.
    So, merge by gunzip and gzip again...
    """
    input:
        fofn = 'input/fastq/{mrg_sample}_1000.fofn'
    output:
        'input/fastq/{mrg_sample}_1000.fastq.gz'
    log:
        'log/input/fastq/{mrg_sample}_1000.merge.log'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    conda:
        '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_medium']
    resources:
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input),
        threads = lambda wildcards: config['num_cpu_low'] // 2
    shell:
        'pigz -p {params.threads} -d -c {params.fastq_parts} | pigz -p {params.threads} > {output} 2> {log}'


def collect_pacbio_bam_input_parts(wildcards, glob_collect=True):
    """
    :param wildcards:
    :return:
    """
    import os
    readset = wildcards.readset
    subfolder = os.path.join('bam', readset)

    if glob_collect:
        import glob
        pattern = os.path.join('input', subfolder, readset + '.part*.pbn.bam')
        bam_parts = glob.glob(pattern)

        if not bam_parts:
            raise RuntimeError('collect_pacbio_bam_input_parts: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        request_path = checkpoints.create_input_data_download_requests.get(subfolder='bam', readset=readset).output[0]
        data_path = os.path.split(request_path)[0]

        checkpoint_wildcards = glob_wildcards(os.path.join(request_path, readset + '.{part_num}.request'))

        bam_parts = expand(
            os.path.join(data_path, readset + '.{part_num}.pbn.bam'),
            part_num=checkpoint_wildcards.part_num
        )

    return bam_parts


rule write_bam_input_parts_fofn:
    input:
        req_dir = 'input/bam/{readset}/requests',
        bam_parts = collect_pacbio_bam_input_parts
    output:
        fofn = 'input/bam/{readset}_1000.pbn.fofn'
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
    """
    switch to pbmerge ; bamtools merge
    seems to create corrupt data blocks
    every now and then
    """
    input:
        fofn = 'input/bam/{mrg_sample}_1000.pbn.fofn'
    output:
        'input/bam/{mrg_sample}_1000.pbn.bam'
    log:
        'log/input/bam/{mrg_sample}_1000.mrg.log'
    benchmark:
        'rsrc/input/bam/{mrg_sample}_1000.mrg.rsrc'
    wildcard_constraints:
        mrg_sample = CONSTRAINT_PARTS_PBN_INPUT_SAMPLES
    conda:
         '../environment/conda/conda_pbtools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 if (attempt <= 1 and '-ccs' in wildcards.mrg_sample) else 24 * attempt
    params:
        bam_parts = lambda wildcards, input: load_fofn_file(input)
    shell:
        'pbmerge {params.bam_parts} > {output} 2> {log}'


rule chs_child_filter_to_100x:
    """
    This one sample has ~200x coverage, and cannot be processed by flye
    Hard-code for now as this is not required for any other input sample
    """
    input:
        'input/bam/HG00514_hgsvc_pbsq2-clr_1000.pbn.bam'
    output:
        'input/bam/HG00514_hgsvc_pbsq2-clr_0526.pbn.bam'
    log:
        'log/input/bam/HG00514_hgsvc_pbsq2-clr_0526.sampling.log'
    benchmark:
        'rsrc/input/bam/HG00514_hgsvc_pbsq2-clr_0526.sampling.rsrc'
    message: 'DEPRECATED: downsampling for sample HG00514 was only done for early chemistry read data'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'bamtools filter -length ">=30000" -in {input} -out {output} &> {log}'


def collect_strandseq_libraries(wildcards, glob_collect=True, caller='snakemake'):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    import os
    source_path = 'input/fastq/{sseq_reads}/{lib_id}.fastq.gz'

    if glob_collect:
        import glob
        pattern = source_path.replace('{lib_id}', '*')
        pattern = pattern.format(**dict(wildcards))
        sseq_fastq = glob.glob(pattern)

        if not sseq_fastq:
            if caller == 'snakemake':
                tmp = dict(wildcards)
                sseq_fastq = []
                sseq_libs, sseq_lib_ids = get_strandseq_library_info(wildcards.sseq_reads)
                for sseq_lib in sseq_libs:
                    tmp['lib_id'] = sseq_lib + '_1'
                    sseq_fastq.append(source_path.format(**tmp))
                    tmp['lib_id'] = sseq_lib + '_2'
                    sseq_fastq.append(source_path.format(**tmp))
            else:
                raise RuntimeError('collect_strandseq_libraries: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        checkpoint_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sseq_reads).output[0]

        glob_pattern = os.path.join(checkpoint_dir, '{lib_id}.request')

        checkpoint_wildcards = glob_wildcards(glob_pattern)

        sseq_fastq = expand(
            source_path,
            sseq_reads=wildcards.sseq_reads,
            lib_id=checkpoint_wildcards.lib_id
            )
    assert sseq_fastq, 'collect_strandseq_libraries >> returned empty output: {}'.format(wildcards)
    return sseq_fastq


rule merge_strandseq_libraries:
    """
    To have a simple way of incorporating the sseq_reads
    wildcard into the workflow, create this file listing
    to be referred to downstream
    """
    input:
        sseq_libs = collect_strandseq_libraries
    output:
        'input/fastq/{sseq_reads}.fofn'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_SAMPLES
    run:
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(input.sseq_libs)))


def collect_short_read_input_parts(wildcards, glob_collect=True):
    import os
    source_path = 'input/fastq/{readset}/{lib_prefix}_{lib_id}_{mate}.fastq.gz'

    lib_prefix = wildcards.readset.rsplit('_', 1)[0]

    if glob_collect:
        import glob
        custom_wildcards = {
            'readset': wildcards.readset,
            'lib_prefix': lib_prefix,
            'mate': wildcards.mate
        }
        pattern = source_path.replace('{lib_id}', '*')
        pattern = pattern.format(**custom_wildcards)
        short_fastq = glob.glob(pattern)

        if not short_fastq:
            raise RuntimeError('collect_short_read_input_parts: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        checkpoint_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.readset).output[0]

        fix_mate_pattern = '_'.join([lib_prefix, '{lib_id}', wildcards.mate])

        glob_pattern = os.path.join(checkpoint_dir, fix_mate_pattern + '.request')

        checkpoint_wildcards = glob_wildcards(glob_pattern)

        short_fastq = expand(
            source_path,
            readset=[wildcards.readset] * len(checkpoint_wildcards.lib_id),
            lib_prefix=[lib_prefix] * len(checkpoint_wildcards.lib_id),
            lib_id=checkpoint_wildcards.lib_id,
            mate=[wildcards.mate] * len(checkpoint_wildcards.lib_id)
        )

    return short_fastq


rule write_short_read_input_fofn:
    input:
        short_reads = collect_short_read_input_parts
    output:
        'input/fastq/{readset}_{mate}.fofn'
    wildcard_constraints:
        readset = CONSTRAINT_SHORT_READ_INPUT_SAMPLES
    run:
        try:
            validate_checkpoint_output(input.short_reads)
            short_fastq = input.short_reads
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            short_fastq = collect_short_read_input_parts(wildcards, glob_collect=True)

        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(short_fastq)))


rule merge_partial_short_read_input_samples:
    input:
         fofn = 'input/fastq/{readset}_{mate}.fofn'
    output:
        'input/fastq/{readset}_{mate}.fastq.gz'
    log:
        'log/input/fastq/{readset}_{mate}.merge.log'
    benchmark:
        'rsrc/input/fastq/{readset}_{mate}.merge.rsrc'
    wildcard_constraints:
        readset = CONSTRAINT_PARTS_SHORT_READ_INPUT_SAMPLES
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input),
        threads = lambda wildcards: config['num_cpu_low'] // 2
    shell:
         'pigz -p {params.threads} -d -c {params.fastq_parts} | pigz -p {threads} > {output} 2> {log}'


rule relink_complete_short_read_input_samples:
    """
    TODO: fix download path for complete samples
    """
    input:
        fofn = 'input/fastq/{readset}_{mate}.fofn'
    output:
        'input/fastq/{readset}_{mate}.fastq.gz'
    log:
        'log/input/fastq/{readset}_{mate}.merge.log'
    benchmark:
        'rsrc/input/fastq/{readset}_{mate}.merge.rsrc'
    wildcard_constraints:
        readset = CONSTRAINT_COMPLETE_SHORT_READ_INPUT_SAMPLES
    params:
        fastq_parts = lambda wildcards, input: load_fofn_file(input)
    shell:
        'ln --symbolic --relative {params.fastq_parts} {output}'

