
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    shell:
        "md5sum {input} > {output}"


rule samtools_index_cram_alignment:
    input:
        cram = '{filepath}.cram'
    output:
        crai = '{filepath}.cram.crai'
    threads: config['num_cpu_low']
    shell:
        "samtools index -@ {threads} {input.cram}"


rule samtools_index_bam_alignment:
    input:
        bam = '{filepath}'
    output:
        bai = '{filepath}.bai'
    threads: config['num_cpu_low']
    shell:
        "samtools index -@ {threads} {input.bam}"


rule pb_bamtools_index_bam_alignment:
    input:
        pbn_bam = '{filepath}.pbn.bam'
    output:
        pbi = '{filepath}.pbn.bam.pbi'
    conda:
        config['conda_env_pbtools']
    shell:
        'pbindex {input}'


rule samtools_convert_sam_to_bam:
    input:
        sam = '{filepath}.sam'
    output:
        bam = '{filepath}.sam.bam'
    wildcard_constraints:
        filepath = '[\w\-\/]+'
    threads: config['num_cpu_low']
    shell:
        "samtools view -o {output.bam} -b -@ {threads} {input.sam}"


rule samtools_position_sort_bam_alignment:
    input:
        unsorted_bam = '{filepath}.sam.bam'
    output:
        sorted_bam = '{filepath}.psort.sam.bam'
    wildcard_constraints:
        filepath = '[\w\-\/]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 20 * 1024,  # 20G per CPU
        mem_total_mb = config['num_cpu_low'] * (20 * 1024)
    params:
        mem_per_thread = '20G'
    run:
        exec = 'samtools sort'
        exec += ' -m {params.mem_per_thread}'
        exec += ' --threads {threads}'
        exec += ' -o {output.sorted_bam}'
        exec += ' {input.unsorted_bam}'
        shell(exec)


rule samtools_index_fasta:
    input:
        '{filepath}.fasta'
    output:
        '{filepath}.fasta.fai'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    shell:
        'bcftools index --tbi {input}'


checkpoint create_assembly_sequence_files:
    input:
        '{folder_path}/{reference}.fasta.fai'
    output:
        directory('{folder_path}/{reference}/sequences')
    run:
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
        '{folder_path}/bwa_index/{reference}.amb',
        '{folder_path}/bwa_index/{reference}.ann',
        '{folder_path}/bwa_index/{reference}.bwt',
        '{folder_path}/bwa_index/{reference}.pac',
        '{folder_path}/bwa_index/{reference}.sa'
    log:
        'log/{folder_path}/bwa_index/{reference}.log'
    benchmark:
        'run/{folder_path}/bwa_index/{reference}.rsrc'
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 6144
    params:
        prefix = lambda wildcards, output: output[0].split('.')[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


def collect_strandseq_alignments(wildcards):
    """
    """
    individual, project, platform_spec = wildcards.sts_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')
    try:
        bioproject = config['strandseq_to_bioproject'][wildcards.sts_reads]
    except KeyError:
        bioproject = config['strandseq_to_bioproject'][wildcards.sts_reads.rsplit('_', 1)[0]]

    requests_dir = checkpoints.create_bioproject_download_requests.get(individual=individual, bioproject=bioproject).output[0]

    search_path = os.path.join(requests_dir, '{individual}_{project}_{platform_spec}_{lib_id}_{run_id}_1.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    if hasattr(wildcards, 'reference'):
        reference = wildcards.reference
    else:
        # assume squashed assembly
        reference = wildcards.sample + '_sqa-' + wildcards.assembler
    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}/{bioproject}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.psort.mdup.sam.bam{ext}',
        reference=wildcards.reference,
        individual=individual,
        bioproject=bioproject,
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
        with open(file_path, 'r') as info:
            seq_len = info.readline().split()[1]
            try:
                _ = int(seq_len)
            except ValueError:
                raise AssertionError('Extracted seq. length is not an integer: {} / {}'.format(seq_len, file_path))
    return seq_len