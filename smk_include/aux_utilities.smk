
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    conda:
         '../environment/conda/conda_shelltools.yml'
    shell:
        "md5sum {input} > {output}"


rule samtools_index_bam_alignment:
    input:
        bam = '{filepath}'
    output:
        bai = '{filepath}.bai'
    benchmark:
        'run/{{filepath}}.idx-bai.t{}.rsrc'.format(config['num_cpu_low'])
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 8 * attempt
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
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 8 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    shell:
        'pbindex {input} &> {log}'


rule pb_bam2x_dump_fastq:
    input:
        pbn_bam = 'input/bam/{folder_path}/{pbn_sample}{sample_type}.pbn.bam',
        pbn_idx = 'input/bam/{folder_path}/{pbn_sample}{sample_type}.pbn.bam.pbi'
    output:
        'input/fastq/{folder_path}/{pbn_sample}{sample_type}.fastq.gz'
    log:
        'log/input/bam/{folder_path}/{pbn_sample}{sample_type}.dump.log'
    benchmark:
        'run/input/bam/{folder_path}/{pbn_sample}{sample_type}.dump.rsrc'
    wildcard_constraints:
        pbn_sample = '(' + '|'.join(config['partial_pbn_samples'] + config['complete_pbn_samples']) + ')',
        sample_type = '(_[0-9x]+|\.part[0-9]+)'
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 if attempt <= 1 else 12 * attempt
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
        pbn_hap_reads = '(' + '|'.join(config['partial_pbn_samples'] + config['complete_pbn_samples']) + ')_[0-9]+',
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
        'input/fastq/complete/{file_name}.fastq.gz'
    output:
        'input/fasta/complete/{file_name}.fasta'
    log:
        'log/input/fastq/complete/{file_name}.fa-dump.log'
    benchmark:
        'run/input/fastq/complete/{file_name}.fa-dump.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 + 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 + 4096 * attempt
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_dir = config['script_dir'],
        buffer_limit = 8  # gigabyte
    shell:
        '{params.script_dir}/dump_shasta_fasta.py --debug --buffer-size {params.buffer_limit} '
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
        script_dir = config['script_dir'],
        buffer_limit = 4  # gigabyte
    shell:
        '{params.script_dir}/dump_shasta_fasta.py --debug --buffer-size {params.buffer_limit} '
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
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule generate_bwa_index:
    input:
        reference = '{folder_path}/{reference}.fasta'
    output:
        '{folder_path}/{reference}/bwa_index/{reference}.amb',
        '{folder_path}/{reference}/bwa_index/{reference}.ann',
        '{folder_path}/{reference}/bwa_index/{reference}.bwt',
        '{folder_path}/{reference}/bwa_index/{reference}.pac',
        '{folder_path}/{reference}/bwa_index/{reference}.sa'
    log:
        'log/{folder_path}/bwa_index/{reference}.log'
    benchmark:
        'run/{folder_path}/bwa_index/{reference}.rsrc'
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 6144,
        runtime_hrs = 3
    conda:
         '../environment/conda/conda_biotools.yml'
    params:
        prefix = lambda wildcards, output: output[0].split('.')[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


rule singularity_pull_container:
    input:
        'output/check_files/environment/singularity_version.ok'
    output:
        'output/container/{hub}/{repo}/{tool}_{version}.sif'
    log:
        'log/output/container/{hub}/{repo}/{tool}_{version}.pull.log'
    #envmodules:
       #    config['env_module_singularity']
    params:
        pull_folder = lambda wildcards: os.path.join(os.getcwd(), 'output', 'container', wildcards.hub, wildcards.repo),
        singularity_module = config['env_module_singularity']
    shell:
        'module load {params.singularity_module} ; '
        'SINGULARITY_PULLFOLDER={params.pull_folder} singularity pull '
            '{wildcards.hub}://{wildcards.repo}/{wildcards.tool}:{wildcards.version} &> {log} ; '
        'module unload {params.singularity_module}'


def collect_strandseq_alignments(wildcards):
    """
    """
    individual, project, platform_spec = wildcards.sts_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')

    requests_dir = checkpoints.create_bioproject_download_requests.get(sts_reads=wildcards.sts_reads).output[0]

    search_path = os.path.join(requests_dir, '{individual}_{project}_{platform_spec}_{lib_id}_{run_id}_1.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.psort.mdup.sam.bam{ext}',
        reference=wildcards.reference,
        individual=individual,
        sts_reads=wildcards.sts_reads,
        project=project,
        platform=platform,
        lib_id=checkpoint_wildcards.lib_id,
        ext=['', '.bai'])
    return sorted(bam_files)


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