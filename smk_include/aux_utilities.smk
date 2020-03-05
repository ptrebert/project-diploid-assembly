
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    conda:
         '../environment/conda/conda_shelltools.yml'
    shell:
        "md5sum {input} > {output}"


rule gunzip_fastq_copy:
    """
    exists only for Shasta input
    """
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.fastq'
    conda:
         '../environment/conda/conda_shelltools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 12 * attempt
    shell:
        "gzip -dc {input} > {output}"


rule samtools_index_bam_alignment:
    input:
        bam = '{filepath}.bam'
    output:
        bai = '{filepath}.bam.bai'
    benchmark:
        'run/{{filepath}}.idx-bai.t{}.rsrc'.format(config['num_cpu_low'])
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 16 * attempt
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        "samtools index -@ {threads} {input.bam}"


rule pb_bamtools_index_bam_alignment:
    input:
        pbn_bam = '{filepath}.pbn.bam'
    output:
        pbi = '{filepath}.pbn.bam.pbi'
    log:
        'log/{filepath}.create-pbi.log'
    benchmark:
        'run/{filepath}.create-pbi.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 16 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    shell:
        'pbindex {input} &> {log}'


rule pb_bam2x_dump_fastq:
    input:
        pbn_bam = 'input/bam/{pbn_sample}_{sampling}.pbn.bam',
        pbn_idx = 'input/bam/{pbn_sample}_{sampling}.pbn.bam.pbi'
    output:
        'input/fastq/{pbn_sample}_{sampling}.fastq.gz'
    log:
        'log/input/bam/{pbn_sample}_{sampling}.dump.log'
    benchmark:
        'run/input/bam/{pbn_sample}_{sampling}.dump.rsrc'
    wildcard_constraints:
        pbn_sample = CONSTRAINT_ALL_PBN_INPUT_SAMPLES,
        sampling = '[0-9x]+'
    resources:
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0]
    shell:
        'bam2fastq -c 5 -o {params.out_prefix} {input.pbn_bam} &> {log}'


rule pb_bam2x_dump_haploid_fastq:
    input:
        pbn_bam = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{pbn_hap_reads}.{hap}.{sequence}.pbn.bam',
        pbn_idx = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{pbn_hap_reads}.{hap}.{sequence}.pbn.bam.pbi',
    output:
        'output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.fastq.gz',
    log:
        'log/output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.dump.log',
    benchmark:
        'run/output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.dump.rsrc',
    wildcard_constraints:
        pbn_hap_reads = CONSTRAINT_ALL_PBN_INPUT_SAMPLES + '_[0-9]+',
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 2 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0]
    shell:
        'bam2fastq -c 5 -o {params.out_prefix} {input.pbn_bam} &> {log}'


rule dump_shasta_fasta:
    """
    Despite the docs saying that (gzipped) fastq can be used
    as input to Shasta, Shasta complains about it...
    Dump to single-line (!) FASTA, only format that is accepted
    (tested with v0.3.0)
    """
    input:
        'input/fastq/{file_name}.fastq.gz'
    output:
        'input/fasta/{file_name}.fasta'
    log:
        'log/input/fastq/{file_name}.fa-dump.log'
    benchmark:
        'run/input/fastq/{file_name}.fa-dump.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 + 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 + 4096 * attempt
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('dump_shasta_fasta.py'),
        buffer_limit = 8  # gigabyte
    shell:
        '{params.script_exec} --debug --buffer-size {params.buffer_limit} '
            ' --input-fq {input} --output-fa {output} &> {log}'


rule dump_shasta_haploid_fasta:
    """
    Despite the docs saying that (gzipped) fastq can be used
    as input to Shasta, Shasta complains about it...
    Dump to single-line (!) FASTA, only format that is accepted
    (tested with v0.3.0)
    """
    input:
        '{file_path}/haploid_fastq/{file_name}.fastq.gz'
    output:
        '{file_path}/haploid_fasta/{file_name}.fasta'
    log:
        'log/{file_path}/haploid_fastq/{file_name}.fa-dump.log'
    benchmark:
        'run/{file_path}/haploid_fastq/{file_name}.fa-dump.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 2 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 + 4096 * attempt
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('dump_shasta_fasta.py'),
        buffer_limit = 4  # gigabyte
    shell:
        '{params.script_exec} --debug --buffer-size {params.buffer_limit} '
            ' --input-fq {input} --output-fa {output} &> {log}'


rule samtools_index_fasta:
    input:
        '{filepath}.fasta'
    output:
        '{filepath}.fasta.fai'
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        'bcftools index --tbi {input}'


checkpoint create_assembly_sequence_files:
    input:
        '{folder_path}/{reference}.fasta.fai'
    output:
        directory('{folder_path}/{reference}/sequences')
    run:
        import os
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                output_path = os.path.join(output_dir, seq_name + '.seq')
                if os.path.isfile(output_path):
                    with open(output_path, 'r') as seq_file:
                        seq_entry = seq_file.read().strip()
                        if seq_entry == line.strip():
                            continue
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule generate_bwa_index:
    input:
        reference = '{folder_path}/{file_name}.fasta'
    output:
        '{folder_path}/{file_name}/bwa_index/{file_name}.amb',
        '{folder_path}/{file_name}/bwa_index/{file_name}.ann',
        '{folder_path}/{file_name}/bwa_index/{file_name}.bwt',
        '{folder_path}/{file_name}/bwa_index/{file_name}.pac',
        '{folder_path}/{file_name}/bwa_index/{file_name}.sa'
    log:
        'log/{folder_path}/bwa_index/{file_name}.log'
    benchmark:
        'run/{folder_path}/bwa_index/{file_name}.rsrc'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    conda:
         '../environment/conda/conda_biotools.yml'
    params:
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


rule singularity_pull_container:
    input:
        'output/check_files/environment/singularity_version.ok'
    output:
        'output/container/{hub}/{repo}/{tool}_{version}.sif'
    log:
        'log/output/container/{hub}/{repo}/{tool}_{version}.pull.log'
    envmodules:
        config['env_module_singularity']
    params:
        pull_folder = lambda wildcards: os.path.join(os.getcwd(), 'output', 'container', wildcards.hub, wildcards.repo),
        singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
    shell:
        '{params.singularity}'
        'SINGULARITY_PULLFOLDER={params.pull_folder} singularity pull '
            '{wildcards.hub}://{wildcards.repo}/{wildcards.tool}:{wildcards.version} &> {log}'


def collect_strandseq_alignments(wildcards, glob_collect=False):
    """
    """
    source_path = 'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup.sam.bam{ext}'

    individual, project, platform_spec = wildcards.sts_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')

    if glob_collect:
        import glob
        source_path = source_path.replace('{ext}', '*')
        source_path = source_path.replace('{lib_id}', '*')
        pattern = source_path.format(**{'reference': wildcards.reference,
                                        'sts_reads': wildcards.sts_reads,
                                        'individual': individual,
                                        'project': project,
                                        'platform': platform,
                                        'spec': spec})
        bam_files = glob.glob(pattern)
        if not bam_files:
            raise RuntimeError('collect_strandseq_alignments: no files collected with pattern {}'.format(pattern))

    else:
        requests_dir = checkpoints.create_input_data_download_requests.get(subfolder='fastq', readset=wildcards.sts_reads).output[0]

        # this is a bit tricky given that there are different varieties of Strand-seq libraries
        glob_pattern = '_'.join([individual, project, platform + '-{spec,[0-9a-z]+}', '{lib_id}'])

        if wildcards.sts_reads in CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES:
            glob_pattern += '_{run_id,[A-Z0-9]+}_1.request'
        else:
            glob_pattern += '_1.request'
        search_path = os.path.join(requests_dir, glob_pattern)

        checkpoint_wildcards = glob_wildcards(search_path)

        bam_files = expand(
            source_path,
            reference=wildcards.reference,
            individual=individual,
            sts_reads=wildcards.sts_reads,
            project=project,
            platform=platform,
            spec=spec,
            lib_id=checkpoint_wildcards.lib_id,
            ext=['', '.bai'])

    return sorted(bam_files)


def validate_checkpoint_output(check_output, expected_type='list_of_files'):
    """
    Because of github issues #55 and #142, this function exists to make
    sure the output returned from a checkpoint evaluation conforms to
    expectation (so far always list of file paths)

    :param check_output:
    :param expected_type:
    :return:
    """
    if expected_type == 'list_of_files':
        if isinstance(check_output, list):
            num_of_items = len(check_output)
            if num_of_items == 0:
                raise ValueError('Checkpoint evaluation resulted in empty list of files')
            for item in check_output:
                if os.path.isdir(item):
                    raise ValueError('Validating checkpoint output w/ {} items '
                                     '- encountered directory: {}'.format(num_of_items, item))
        elif isinstance(check_output, str):
            if os.path.isdir(check_output):
                # this is the reported case, handle in particular
                raise RuntimeError('Checkpoint evaluation failed, '
                                   'returned single str / directory: {}'.format(check_output))
            else:
                raise ValueError('Single string returned from checkpoint evaluation, '
                                 'but is (not yet) a directory path: {}'.format(check_output))
        else:
            raise ValueError('Non-list type received for checkpoint '
                             'output validation: {} / {}'.format(type(check_output), str(check_output)))
    else:
        raise NotImplementedError('aux_utilities::validate_checkpoint_output called with '
                                  'unsupported expected type: {}'.format(expected_type))
    return


def load_preset_file(wildcards, input):
    """
    Save load for parameters from files that do not
    exist during a dry run, leading to a FileNotFoundError
    raised by Snakemake. Replaces construct

    params:
        preset = lambda wildcards, input: open(input.preset).read().strip()

    with

    params:
        preset = load_preset_file

    """
    if not hasattr(input, 'preset'):
        raise AttributeError('Input does not have a "preset" attribute: {}'.format(input))
    file_path = input.preset
    if not os.path.isfile(file_path):
        preset = 'PRESET-DRY-RUN'
    else:
        with open(file_path, 'r') as dump:
            preset = dump.read().strip()
            assert preset, 'Empty preset file: {}'.format(file_path)
    return preset


def load_fofn_file(input, prefix='', sep=' '):
    """
    Save load for list of filenames from a file that
    does not exist during a dry run, leading to a
    FileNotFoundError raised by Snakemake.

    Replaces construct

    params:
        preset = lambda wildcards, input: open(input.fofn).read().strip()

    with

    params:
        preset = lambda wildcards, input: load_fofn_file(input, prefix_string, separator_string)
    """
    if not hasattr(input, 'fofn'):
        raise AttributeError('Input does not have "fofn" attribute: {}'.format(input))
    file_path = input.fofn
    if not os.path.isfile(file_path):
        file_list = ['FOFN-DRY-RUN']
    else:
        with open(file_path, 'r') as dump:
            file_list = sorted([l.strip() for l in dump.readlines()])
            assert file_list, 'Empty fofn file: {}'.format(file_path)
    file_list = prefix + sep.join(file_list)
    return file_list


def load_seq_length_file(wildcards, input):
    """
    Follows same idea as two function above
    """
    if not hasattr(input, 'seq_info'):
        raise AttributeError('Input does not have "seq_info" attribute: {}'.format(input))
    file_path = input.seq_info
    if not os.path.isfile(file_path):
        seq_len = 'SEQLEN-DRY-RUN'
    else:
        seq_len = 0
        with open(file_path, 'r') as seq_info:
            for line in seq_info:
                if line.startswith('#'):
                    continue
                length = line.split()[1]
                try:
                    length = int(length)
                    seq_len += length
                except ValueError:
                    raise ValueError('Extracted seq. length is not an integer: {} / {}'.format(line.strip(), file_path))
    return seq_len


def find_script_path(script_name, subfolder=''):
    """
    Find full path to script to be executed. Function exists
    to avoid config parameter "script_dir"

    :param script_name:
    :param subfolder:
    :return:
    """
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path
